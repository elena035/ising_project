#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include "utils.h"

// ============================================================================
// SEZIONE 1: INIZIALIZZAZIONE E LETTURA DATI
// ============================================================================

void leggi_parametri(const char *filename, parametri *p) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "ERRORE: Impossibile aprire il file di configurazione %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // 1. Lettura costanti fisiche
    if (fscanf(file, "%*s %lf", &p->k_B) != 1) {
        fprintf(stderr, "ERRORE: impossibile leggere k_B dal file di configurazione.\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(file, "%*s %lf", &p->h) != 1) {
        fprintf(stderr, "ERRORE: impossibile leggere h dal file di configurazione.\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(file, "%*s %lf", &p->J) != 1) {
        fprintf(stderr, "ERRORE: impossibile leggere J dal file di configurazione.\n");
        exit(EXIT_FAILURE);
    }

    // 2. Lettura num_betas e array dei beta
    if (fscanf(file, "%*s %d", &p->num_betas) != 1) {
        fprintf(stderr, "ERRORE: impossibile leggere num_betas dal file di configurazione.\n");
        exit(EXIT_FAILURE);
    }
    
    p->beta_values = malloc(p->num_betas * sizeof(double));
    
    if (fscanf(file, "%*s") == EOF) { // Salta l'etichetta "betas:"
        fprintf(stderr, "ERRORE: fine inaspettata del file prima dei valori beta.\n");
        exit(EXIT_FAILURE);
    } 
    
    for (int i = 0; i < p->num_betas; i++) {
        if (fscanf(file, "%lf", &p->beta_values[i]) != 1) {
            fprintf(stderr, "ERRORE: impossibile leggere il valore beta all'indice %d.\n", i);
            exit(EXIT_FAILURE);
        }
    }

    // 3. Lettura parametri operativi e geometrici
    if (fscanf(file, "%*s %d", &p->frame_freq) != 1) {
        fprintf(stderr, "ERRORE: impossibile leggere frame_freq dal file di configurazione.\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(file, "%*s %d", &p->meas_freq) != 1) {
        fprintf(stderr, "ERRORE: impossibile leggere meas_freq dal file di configurazione.\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(file, "%*s %d", &p->meas_sweeps) != 1) {
        fprintf(stderr, "ERRORE: impossibile leggere meas_sweeps dal file di configurazione.\n");
        exit(EXIT_FAILURE);
    }
    if (fscanf(file, "%*s %d", &p->eq_sweeps) != 1) {
        fprintf(stderr, "ERRORE: impossibile leggere eq_sweeps dal file di configurazione.\n");
        exit(EXIT_FAILURE);
    }

    // 4. Lettura num_L e array degli L
    if (fscanf(file, "%*s %d", &p->num_L) != 1) {
        fprintf(stderr, "ERRORE: impossibile leggere num_L dal file di configurazione.\n");
        exit(EXIT_FAILURE);
    }
    
    
    p->L_values = malloc(p->num_L * sizeof(int));

    if (fscanf(file, "%*s") == EOF) { // Salta l'etichetta "L_values"
        fprintf(stderr, "ERRORE: fine inaspettata del file prima dei valori L.\n");
        exit(EXIT_FAILURE);
    } 
    
    for (int i = 0; i < p->num_L; i++) {
        if (fscanf(file, "%d", &p->L_values[i]) != 1) {
            fprintf(stderr, "ERRORE: impossibile leggere il valore L all'indice %d.\n", i);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
    printf("Parametri caricati correttamente da: %s\n", filename);
}


void init_rng(int rank) {
    srand(time(NULL) + rank);
}


// ============================================================================
// SEZIONE 2: GESTIONE INPUT/OUTPUT E SALVATAGGIO
// ============================================================================



void inizializza_file_csv(const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "ERRORE: Impossibile creare il file dati %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fprintf(file, "L,beta,<M>,<M2>,errore su <M>\n");
    fclose(file);

    printf("File dati %s inizializzato.\n", filename);
}



void salva_misura_csv(const char *filename, int L, double beta, double mag_media, double mag2_media, double errore_mag) {
    FILE *file = fopen(filename, "a");
    if (file == NULL) {
        fprintf(stderr, "ERRORE: Impossibile aprire il file dati %s per il salvataggio\n", filename);
        exit(EXIT_FAILURE);
    }

    fprintf(file, "%d,%f,%f,%f,%f\n", L, beta, mag_media, mag2_media, errore_mag);
    fclose(file);
}



void salva_configurazione(FILE *file, int *reticolo, int L, double beta, int sweep) {
    fprintf(file, "L: %d BETA: %.3f SWEEP: %d\n", L, beta, sweep);
    for (int j = 0; j < L; j++) {
        for (int i = 0; i < L; i++) {
            int k = i + j * L;
            fprintf(file, "%2d ", reticolo[k]); 
        }
        fprintf(file, "\n");
    }
}


// ============================================================================
// SEZIONE 3: COMUNICAZIONE MPI
// ============================================================================


void raccordo_dominio(int *reticolo, int L, int L_local, int rank, int size) {
    // identificazione topologia MPI
    int vicino_sopra = (rank - 1 + size) % size;
    int vicino_sotto = (rank + 1) % size;

    // Impostazione dei puntatori alle righe di scambio
    int *riga_fantasma_top = &reticolo[0 * L];                
    int *prima_riga_reale  = &reticolo[1 * L];                
    int *ultima_riga_reale = &reticolo[L_local * L];          
    int *riga_fantasma_bot = &reticolo[(L_local + 1) * L];    

    // Scorrimento verso l'alto: invio la prima reale in alto, ricevo nel bordo fantasma inferiore
    MPI_Sendrecv(prima_riga_reale, L, MPI_INT, vicino_sopra, 0,
                 riga_fantasma_bot, L, MPI_INT, vicino_sotto, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Scorrimento verso il basso: invio l'ultima riga reale in basso, ricevo nel bordo fantasma superiore
    MPI_Sendrecv(ultima_riga_reale, L, MPI_INT, vicino_sotto, 1,
                 riga_fantasma_top, L, MPI_INT, vicino_sopra, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
