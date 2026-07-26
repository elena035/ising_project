#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "utils.h"
#include "lattice.h"
#include "metropolis.h"
#include "measurements.h"
#include "validazione.h"

int main(int argc, char *argv[]) {
    int rank, size;
    
    // --- 1. INIZIALIZZAZIONE MPI ---
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3) {
        if (rank == 0) {
            printf("ERRORE: Uso: mpirun -np <N> %s <input.txt> <output.csv>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    const char *file_input = argv[1]; 
    const char *file_output = argv[2]; 
    
    parametri p; 
    leggi_parametri(file_input, &p);

    // --- CONTROLLO ROBUSTEZZA: DIVISIONE ESATTA DEL RETICOLO ---
    // Controlliamo subito che tutti i valori di L siano divisibili per il numero di processi
    for (int i = 0; i < p.num_L; i++) {
        if (p.L_values[i] % size != 0) {
            if (rank == 0) {
                printf("ERRORE CRITICO: Il reticolo L = %d non è divisibile esattamente per %d processi.\n", p.L_values[i], size);
            }
            MPI_Finalize();
            return 1;
        }
    }

    // --- STAMPA INTESTAZIONE E PARAMETRI (Solo Rank 0) ---
    if (rank == 0) {
        printf("\n================================================================================\n");
        printf("           SIMULAZIONE MODELLO DI ISING 2D (Metropolis MPI)\n");
        printf("================================================================================\n");
        printf("Questo programma simula l'evoluzione del modello di Ising su reticoli 2D per\n");
        printf("studiare la transizione di fase. Per ogni combinazione di (L, beta), il codice\n");
        printf("calcola la magnetizzazione media e il suo errore standard (std).\n\n");
        
        printf("Parametri caricati da: %s\n", file_input);
        
        printf(" - Taglie reticolo (L): [");
        for (int i = 0; i < p.num_L; i++) {
            printf("%d%s", p.L_values[i], (i < p.num_L - 1) ? ", " : "]\n");
        }
        
        printf(" - Temperature (beta):  [");
        for (int j = 0; j < p.num_betas; j++) {
            printf("%.3f%s", p.beta_values[j], (j < p.num_betas - 1) ? ", " : "]\n");
        }
        printf("================================================================================\n\n");
    }

    if (rank == 0) esegui_test_validazione(&p); 
    
    init_rng(rank); 

    // --- PULIZIA E PREPARAZIONE CARTELLE GLOBALI ---
    if (rank == 0) {
        inizializza_file_csv(file_output); 
        mkdir("simulazioni", 0777);
    }

    // --- 2. CICLO SULLE TAGLIE (L) ---
    for(int i = 0; i < p.num_L; i++) {
        
        // Timer per i tempi di calcolo e di comunicazione
        double tempo_calcolo_totale = 0.0;
        double tempo_comunicazione_totale = 0.0;

        p.L = p.L_values[i]; 
        int L = p.L;
        int L_local = L / size; 

        if (rank == 0) {
            printf("\n=== Inizio simulazione per L = %d ===\n", L);
            char dir_L[100];
            sprintf(dir_L, "simulazioni/L_%d", L);
            mkdir(dir_L, 0777);
        }

        int *reticolo = costruisci_reticolo(&p, L_local);
        vicini *lista_vicini = definisci_vicini(&p, L_local);

        int *indici_rossi = NULL;
        int *indici_neri = NULL;
        costruisci_liste_scacchiera(L, L_local, rank, &indici_rossi, &indici_neri);

        raccordo_dominio(reticolo, L, L_local, rank, size);

        // --- 3. CICLO SULLE TEMPERATURE (BETA) ---
        for(int j = 0; j < p.num_betas; j++) {
            double beta = p.beta_values[j];
            double probabilities[9]; 
            init_probabilities(beta, p.J, probabilities);

            FILE *file_evoluzione = NULL;
            char nome_traiettoria[256]; 
            
            if (rank == 0) {
                char nome_evoluzione[256];
                sprintf(nome_evoluzione, "simulazioni/L_%d/evoluzione_L_%d_beta_%.3f.csv", L, L, beta);
                file_evoluzione = fopen(nome_evoluzione, "w");
                if (file_evoluzione != NULL) {
                    fprintf(file_evoluzione, "sweep,magnetizzazione_su_spin\n");
                }
                
                sprintf(nome_traiettoria, "simulazioni/L_%d/traiettoria_L_%d_beta_%.3f.txt", L, L, beta);
                FILE *ftraj = fopen(nome_traiettoria, "w");
                if (ftraj != NULL) fclose(ftraj);
            }

            // --- 4. TERMALIZZAZIONE ---
            for (int eq = 0; eq < p.eq_sweeps; eq++) {
                metropolis_sweep_completo(reticolo, lista_vicini, probabilities, indici_rossi, indici_neri, L, L_local, rank, size, &tempo_calcolo_totale, &tempo_comunicazione_totale);
            }

            double sum_m = 0.0, sum_m2 = 0.0;
            int campioni = 0;

            int *reticolo_globale = NULL;
            if (rank == 0) {
                reticolo_globale = malloc(L * L * sizeof(int));
            }

            int max_campioni = (p.meas_sweeps / p.meas_freq) + 1;
            int *buffer_meas = NULL;
            double *buffer_m_ist = NULL;
            int buffer_idx = 0;
            
            if (rank == 0 && file_evoluzione != NULL) {
                buffer_meas = malloc(max_campioni * sizeof(int));
                buffer_m_ist = malloc(max_campioni * sizeof(double));
            }

            // --- 5. MISURAZIONE ---
            for (int meas = 1; meas <= p.meas_sweeps; meas++) {
                
                metropolis_sweep_completo(reticolo, lista_vicini, probabilities, indici_rossi, indici_neri, L, L_local, rank, size, &tempo_calcolo_totale, &tempo_comunicazione_totale);

                // CAMPIONAMENTO FISICO
                if (meas % p.meas_freq == 0) {
                    double mag_locale = calcola_magnetizzazione_locale(reticolo, L, L_local);
                    double mag_globale = 0.0;
                    
                  
                    double t_comm_start = MPI_Wtime();
                    MPI_Reduce(&mag_locale, &mag_globale, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                    tempo_comunicazione_totale += (MPI_Wtime() - t_comm_start);

                    if (rank == 0) {
                        double m_ist = mag_globale / (L * L);
                        double m_abs = fabs(m_ist);
                        
                        sum_m += m_abs;
                        sum_m2 += m_abs * m_abs;
                        campioni++;
                        
                        if (file_evoluzione != NULL && buffer_meas != NULL && buffer_m_ist != NULL) {
                            if (buffer_idx < max_campioni) {
                                buffer_meas[buffer_idx] = meas;
                                buffer_m_ist[buffer_idx] = m_ist;
                                buffer_idx++;
                            }
                        }
                    }
                }

                // SALVATAGGIO TRAIETTORIE
                if (meas % p.frame_freq == 0) {
                    int elementi_reali_locali = L * L_local;
                    
                    double t_comm_start = MPI_Wtime();
                    MPI_Gather(&reticolo[L], elementi_reali_locali, MPI_INT, reticolo_globale, elementi_reali_locali, MPI_INT, 0, MPI_COMM_WORLD);
                    tempo_comunicazione_totale += (MPI_Wtime() - t_comm_start);
                    
                    if (rank == 0) {
                        FILE *fe = fopen(nome_traiettoria, "a");
                        if (fe != NULL) {
                            salva_configurazione(fe, reticolo_globale, L, beta, meas);
                            fclose(fe);
                        }
                    }
                }
            }

            // CHIUSURA FILE EVOLUZIONE E SALVATAGGIO MEDIE
            if (rank == 0) {
                if (file_evoluzione != NULL) {
                    for (int k = 0; k < buffer_idx; k++) {
                        fprintf(file_evoluzione, "%d,%.6f\n", buffer_meas[k], buffer_m_ist[k]);
                    }
                    fclose(file_evoluzione);
                }
                
                if (buffer_meas != NULL) free(buffer_meas);
                if (buffer_m_ist != NULL) free(buffer_m_ist);
                if (reticolo_globale != NULL) free(reticolo_globale);

                double m_media = sum_m / campioni;
                double m2_media = sum_m2 / campioni;
                double varianza = (m2_media - (m_media * m_media) < 0) ? 0 : (m2_media - (m_media * m_media));
                double std = sqrt(varianza / campioni);
                
                salva_misura_csv(file_output, L, beta, m_media, m2_media, std);
                printf("Beta: %.3f | <|M|>: %.4f ± %.4f\n", beta, m_media, std);
            }
        }

        // SALVATAGGIO TEMPI PER L SPECIFICO
        if (rank == 0) {
            char nome_file[50];
            sprintf(nome_file, "tempi_speedup_%d.csv", size);
            
            FILE *file_tempi = fopen(nome_file, "a");
            if (file_tempi != NULL) {
                fseek(file_tempi, 0, SEEK_END);
                if (ftell(file_tempi) == 0) {
                    
                    fprintf(file_tempi, "L,comp_time,comm_time\n");
                }
                fprintf(file_tempi, "%d,%f,%f\n", L, tempo_calcolo_totale, tempo_comunicazione_totale);
                fclose(file_tempi);
            }
        }

        free(reticolo);
        free(lista_vicini);
        free(indici_rossi);
        free(indici_neri);

    }

    // --- Tutti i processi liberano i vettori ---
    free(p.L_values);
    free(p.beta_values);

    if (rank == 0) {
        printf("\nSimulazione completata con successo.\n");
    }

    MPI_Finalize();
    return 0;
}
