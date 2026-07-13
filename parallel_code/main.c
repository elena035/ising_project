// ============================================================================
// SIMULAZIONE MODELLO DI ISING 2D - ENTRY POINT (MPI)
// ============================================================================

/*
 * Modulo: main.c
 * --------------
 * Entry point del programma. Orchestra il ciclo di vita dell'ambiente MPI,
 * l'inizializzazione delle strutture dati, l'esecuzione dei loop termodinamici
 * (per taglia L e temperatura beta) e la raccolta globale delle osservabili.
 * Gestisce inoltre il profiling dei tempi di esecuzione per l'analisi dello speedup.
 */

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
    
    // ========================================================================
    // FASE 1: INIZIALIZZAZIONE AMBIENTE MPI E SETUP PARAMETRI
    // ========================================================================
    
    // Setup del comunicatore globale
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Validazione degli argomenti a riga di comando
    if (argc != 3) {
        if (rank == 0) {
            printf("\n[ERRORE] Sintassi di avvio non valida.\n");
            printf("Uso previsto: mpirun -np <N> %s <input.txt> <output.csv>\n\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    const char *file_input = argv[1]; 
    const char *file_output = argv[2]; 
    
    // Lettura centralizzata dei parametri fisici e operativi
    parametri p; 
    leggi_parametri(file_input, &p);

    // Stampa dell'intestazione e del recap dei parametri di simulazione (eseguita solo dal Root)
    if (rank == 0) {
        printf("\n================================================================================\n");
        printf("               SIMULAZIONE MODELLO DI ISING 2D (Metropolis MPI)\n");
        printf("================================================================================\n");
        printf(" Questo programma simula il Modello di Ising 2D su un reticolo quadrato di lato L.\n");
	printf(" La simulazione viene effettuata per diverse taglie L e per diversi valori di temperatura inversa beta.\n");
	printf(" L'obiettivo è caratterizzare la transizione di fase che si osserva in prossimità di beta = 0.4406,\n");
	printf(" Misurando la magnetizzazione istantanea per ogni valore di beta utilizzato.\n\n");
        printf(" Algoritmo usato: Metropolis con Domain Decomposition (Red-Black updating).\n\n");
        
        printf(" [CONFIGURAZIONE CARICATA]: %s\n", file_input);
        
        printf(" - Taglie reticolo (L) : [");
        for (int i = 0; i < p.num_L; i++) {
            printf("%d%s", p.L_values[i], (i < p.num_L - 1) ? ", " : "]\n");
        }
        
        printf(" - Temperature (beta)  : [");
        for (int j = 0; j < p.num_betas; j++) {
            printf("%.3f%s", p.beta_values[j], (j < p.num_betas - 1) ? ", " : "]\n");
        }
        printf("================================================================================\n");
    }

    // Esecuzione automatizzata dei test unitari per validare geometria e logica termodinamica
    if (rank == 0) esegui_test_validazione(&p); 
    
    // Inizializzazione del seed pseudo-casuale (scorrelato per ciascun rank MPI)
    init_rng(rank); 

    // Preparazione della gerarchia delle directory per l'esportazione sicura dei log
    if (rank == 0) {
        inizializza_file_csv(file_output); 
        mkdir("simulazioni", 0777);
    }

    // ========================================================================
    // FASE 2: ITERAZIONE SULLE TAGLIE DEL RETICOLO
    // ========================================================================
    
    for(int i = 0; i < p.num_L; i++) {
        
        // Setup dei timer per il profiling dell'efficienza algoritmica
        double tempo_calcolo_totale = 0.0;
        double tempo_comunicazione_totale = 0.0;
        double tempo_totale = 0.0; 
        
	double start;
	double start_simulation = MPI_Wtime();

        p.L = p.L_values[i]; 
        int L = p.L;
        // La taglia locale rappresenta il numero di righe "reali" assegnate a questo processo
        int L_local = L / size; 

        // Generazione delle cartelle di output specifiche per la taglia corrente
        if (rank == 0) {
            char dir_L[100];
            sprintf(dir_L, "simulazioni/L_%d", L);
            mkdir(dir_L, 0777);
            
            printf("\nSimulazione per taglia globale L = %d\n", L);
            printf("---------------------------------------------------------------------\n");
            printf(" %-10s | %-15s | %-15s | %-15s \n", "Beta", "Mag. <|M|>", "std dev");
            printf("---------------------------------------------------------------------\n");
        }


	start = MPI_Wtime();
        // Allocazione dinamica del dominio locale e pre-calcolo delle adiacenze topologiche
        int *reticolo = costruisci_reticolo(&p, L_local);
        vicini *lista_vicini = definisci_vicini(&p, L_local);

        // Classificazione dei siti per l'aggiornamento parallelo indipendente
        int *indici_rossi = NULL;
        int *indici_neri = NULL;
        costruisci_liste_scacchiera(L, L_local, rank, &indici_rossi, &indici_neri);
        
        tempo_calcolo_totale += (MPI_Wtime() - start);

	start = MPI_Wtime();
        // Risoluzione a tempo zero delle ghost cells per preparare i bordi del dominio
        raccordo_dominio(reticolo, L, L_local, rank, size);
        tempo_comunicazione_totale += (MPI_Wtime() - start);

        // ========================================================================
        // FASE 3: ITERAZIONE TERMODINAMICA SULLE TEMPERATURE
        // ========================================================================
        
        for(int j = 0; j < p.num_betas; j++) {
	    start = MPI_Wtime();
            double beta = p.beta_values[j];
            
            // Pre-calcolo dei pesi esponenziali di Boltzmann per ottimizzare il ciclo critico
            double probabilities[9]; 
            init_probabilities(beta, p.J, probabilities);
            tempo_calcolo_totale += (MPI_Wtime() - start);

            FILE *file_evoluzione = NULL;
            char nome_traiettoria[256]; 
            
            // Setup flussi I/O
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

            // --- 3.A: FASE DI TERMALIZZAZIONE ---
            // Decadimento delle condizioni iniziali randomiche per raggiungere l'equilibrio di Boltzmann
	    
            for (int eq = 0; eq < p.eq_sweeps; eq++) {
                metropolis_sweep_completo(reticolo, lista_vicini, probabilities, indici_rossi, indici_neri, 
                                          L, L_local, rank, size, &tempo_calcolo_totale, &tempo_comunicazione_totale);
            }
            

            double sum_m = 0.0, sum_m2 = 0.0;
            int campioni = 0;

            // --- 3.B: FASE DI MISURA E CAMPIONAMENTO ---
            for (int meas = 1; meas <= p.meas_sweeps; meas++) {
                
                // Evoluzione stocastica del sistema
                metropolis_sweep_completo(reticolo, lista_vicini, probabilities, indici_rossi, indici_neri, 
                                          L, L_local, rank, size, &tempo_calcolo_totale, &tempo_comunicazione_totale);
                
		
                // Estrazione delle osservabili macroscopiche a intervalli prefissati (per minimizzare l'autocorrelazione)
		if (meas % p.meas_freq == 0) {
		    start = MPI_Wtime();
                    double mag_locale = calcola_magnetizzazione_locale(reticolo, L, L_local);
                    double mag_globale = 0.0;
	            tempo_calcolo_totale += (MPI_Wtime() - start);
                    
                    // Riduzione globale della magnetizzazione. Mappato nel timer I/O in quanto latenza asincrona.
                    start = MPI_Wtime();
                    MPI_Reduce(&mag_locale, &mag_globale, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                    tempo_comunicazione_totale += (MPI_Wtime() - start);
                    

                    if (rank == 0) {
                        double m_ist = mag_globale / (L * L);
                        double m_abs = fabs(m_ist);
                        
                        sum_m += m_abs;
                        sum_m2 += m_abs * m_abs;
                        campioni++;
                        
                        // Log dinamico dell'evoluzione della magnetizzazione
                        if (file_evoluzione != NULL) { 
                            fprintf(file_evoluzione, "%d,%.6f\n", meas, m_ist); 
                        }
                    }
                }

                // Esportazione della topologia spaziale completa per analisi visive (GIF)
                if (meas % p.frame_freq == 0) {
                    int elementi_reali_locali = L * L_local;
                    int *reticolo_globale = NULL;
                    
                    if (rank == 0) {
                        reticolo_globale = malloc(L * L * sizeof(int));
                    }
                    
                    // Raccolta asimmetrica (Gather) della griglia globale escludendo la prima riga ghost (&reticolo[L])
                    start = MPI_Wtime();
                    MPI_Gather(&reticolo[L], elementi_reali_locali, MPI_INT, reticolo_globale, elementi_reali_locali, MPI_INT, 0, MPI_COMM_WORLD);
                    tempo_comunicazione_totale += (MPI_Wtime() - start);
                    
                    if (rank == 0) {
                        FILE *fe = fopen(nome_traiettoria, "a");
                        if (fe != NULL) {
                            salva_configurazione(fe, reticolo_globale, L, beta, meas);
                            fclose(fe);
                        } 
                        free(reticolo_globale);
                    }
                }


		MPI_Barrier(MPI_COMM_WORLD);
            }

            // Chiusura dei log e calcolo delle medie termodinamiche di ensemble (Solo Root)
            if (rank == 0) {
                 
                if (file_evoluzione != NULL) fclose(file_evoluzione);

                double m_media = sum_m / campioni;
                double m2_media = sum_m2 / campioni;
                
                // Troncamento numerico cautelare per arginare minime imprecisioni hardware in virgola mobile
                double varianza = (m2_media - (m_media * m_media) < 0) ? 0 : (m2_media - (m_media * m_media));
                double std = sqrt(varianza / campioni);

		
                salva_misura_csv(file_output, L, beta, m_media, m2_media, std);
                 
                
                // Stampa in tempo reale dei risultati allineati nella tabella a terminale
                printf(" %-10.3f | %-15.4f | %-15.4f \n", beta, m_media, std);
            }

	    // barriera per sincronizzare i rank prima del prossimo beta
            MPI_Barrier(MPI_COMM_WORLD);
        }

        if (rank == 0) {
            printf("--------------------------------------------------------------------------------\n");
        }

        // ========================================================================
        // FASE 4: REGISTRAZIONE PROFILING PRESTAZIONALE E PULIZIA MEMORIA LOCALE
        // ========================================================================
        
        // 1. Calcolo del tempo totale locale (wall-clock time) per ciascun rank
        double tempo_totale_locale = MPI_Wtime() - start_simulation;

        // 2. Impacchettamento vettoriale dei tempi locali 
        //    Indici: 0 = calcolo, 1 = comunicazione, 2 = totale
        double tempi_locali[3] = {tempo_calcolo_totale, tempo_comunicazione_totale, tempo_totale_locale};
        double tempi_max[3];

        // 3. Riduzione globale: estrazione del tempo massimo (collo di bottiglia) tra tutti i processi
        MPI_Reduce(tempi_locali, tempi_max, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

        // 4. Salvataggio dei tempi di profilazione (eseguito solo dal Root)
        if (rank == 0) {
            char nome_file[50];
            sprintf(nome_file, "tempi_speedup_%d.csv", size);
            
            FILE *file_tempi = fopen(nome_file, "a");
            if (file_tempi != NULL) {
                // Scrittura dell'header se il file è appena stato creato
                fseek(file_tempi, 0, SEEK_END);
                if (ftell(file_tempi) == 0) {
                    fprintf(file_tempi, "L,comp_time,comm_time,total_time\n");
                }
                
                // Scrittura dei tempi massimi globali per l'analisi di speedup
                fprintf(file_tempi, "%d,%f,%f,%f\n", L, tempi_max[0], tempi_max[1], tempi_max[2]);
                
                fclose(file_tempi);
            }
        }

        // Pulizia della memoria deallocando le strutture locali
        free(reticolo);
        free(lista_vicini);
        free(indici_rossi);
        free(indici_neri);
    }

    // ========================================================================
    // FASE 5: CHIUSURA DELL'AMBIENTE MPI
    // ========================================================================
    
    if (rank == 0) {
        free(p.L_values);
        free(p.beta_values);
        printf("\n[SUCCESS] Tutte le iterazioni completate con successo. Uscita sicura.\n");
        printf("================================================================================\n\n");
    }

    MPI_Finalize();
    return 0;
}
