#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "metropolis.h"



// ============================================================================
// SEZIONE 1: FISICA DEL MODELLO E PROBABILITA'
// ============================================================================

void init_probabilities(double beta, double J, double *probabilities) {
        /*
    	Utilizziamo un array di dimensione 9 per mappare direttamente il valore discreto
    	di Delta(E) (che può essere solo 4 o 8 per le mosse svantaggiose) all'indice 
    	dell'array. Questo evita la necessità di istruzioni ed operazioni aggiuntive
    	durante l'esecuzione di Metropolis.
    	*/
        probabilities[4] = exp(-beta * 4.0 * J);
        probabilities[8] = exp(-beta * 8.0 * J);
}

int calcola_delta_E(int *reticolo, vicini *lista_vicini, int k) {
        int somma_vicini = reticolo[lista_vicini[k].up] +
                       reticolo[lista_vicini[k].dw] +
                       reticolo[lista_vicini[k].dx] +
                       reticolo[lista_vicini[k].sx];

    return 2 * reticolo[k] * somma_vicini;
}



// ============================================================================
// SEZIONE 2: DINAMICA E AGGIORNAMENTO
// ============================================================================


void metropolis_mezzo_sweep(int *reticolo, vicini *lista_vicini, double *probabilities, int *indici_colore, int num_indici) {
    
    // Iterazione lineare deterministica per ottimizzare i memory access pattern
    for (int step = 0; step < num_indici; step++) {
        
        // Accesso contiguo all'array degli indici pre-calcolati
        int k = indici_colore[step];

        // Valutazione della variazione di energia locale
        int delta_E = calcola_delta_E(reticolo, lista_vicini, k);

        // Criterio di accettazione di Metropolis
        if (delta_E <= 0) {
            // Transizione energeticamente favorevole o neutra: accettazione garantita
            reticolo[k] = -reticolo[k];
        } else {
            // Transizione sfavorevole: accettazione probabilistica (fattore di Boltzmann)
            // L'estrazione del numero pseudo-casuale è isolata nel branch condizionale 
            // per minimizzare l'overhead computazionale della funzione rand()
            double rnd = (double)rand() / RAND_MAX;
            if (probabilities[delta_E] > rnd) {
                reticolo[k] = -reticolo[k];
            }
        }
    }
}


void metropolis_sweep_completo(int *reticolo, vicini *lista_vicini, double *probabilities, int *indici_rossi, int *indici_neri, int L, int L_local, int rank, int size, double *t_comp, double *t_comm) {
    
    int num_indici = (L * L_local) / 2;
    double inizio, fine;

    // FASE 1: Aggiornamento degli spin rossi (Calcolo puro)
    inizio = MPI_Wtime();
    metropolis_mezzo_sweep(reticolo, lista_vicini, probabilities, indici_rossi, num_indici);
    fine = MPI_Wtime();
    *t_comp += (fine - inizio);

    // FASE 2: Aggiornamento delle ghost cells per i bordi rossi (Comunicazione)
    inizio = MPI_Wtime();
    raccordo_dominio(reticolo, L, L_local, rank, size);
    fine = MPI_Wtime();
    *t_comm += (fine - inizio);

    // FASE 3: Aggiornamento degli spin neri (Calcolo puro dipendente dai rossi)
    inizio = MPI_Wtime();
    metropolis_mezzo_sweep(reticolo, lista_vicini, probabilities, indici_neri, num_indici);
    fine = MPI_Wtime();
    *t_comp += (fine - inizio);

    // FASE 4: Aggiornamento delle ghost cells per i bordi neri (Comunicazione)
    inizio = MPI_Wtime();
    raccordo_dominio(reticolo, L, L_local, rank, size);
    fine = MPI_Wtime();
    *t_comm += (fine - inizio);
}








