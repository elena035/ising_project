#include <stdlib.h>
#include "lattice.h"
#include "utils.h"


// ============================================================================
// SEZIONE 1: GEOMETRIA E MAPPING
// ============================================================================


int calcola_indice_lessicografico(int L, int i, int j) {
    // i scorre sulle colonne, j scorre sulle righe
    int k = i + j * L;
    return k;
}

coordinate calcola_coordinate_cartesiane(int k, int L) {
    coordinate c;
    c.x = k % L; // Resto: posizione nella riga corrente (i)
    c.y = k / L; // Divisione intera: indice della riga (j)
    return c;
}



// ============================================================================
// SEZIONE 2: INIZIALIZZAZIONE RETICOLO
// ============================================================================


int* costruisci_reticolo(parametri *p, int L_local) {
    // Allocazione: L_local righe reali + 2 righe fantasma (ghost cells per MPI)
    int righe_totali = L_local + 2;
    int *reticolo = malloc(righe_totali * p->L * sizeof(int));
    
    // Inizializzazione solo delle righe reali (da j=1 a j=L_local)
    // Le righe j=0 e j=L_local+1 verranno sovrascritte al primo scambio MPI
    for (int j = 1; j <= L_local; j++) {
        for (int i = 0; i < p->L; i++) {
            int k = calcola_indice_lessicografico(p->L, i, j);
            // Assegnazione casuale dello spin: -1 o +1
            reticolo[k] = (rand() % 2) * 2 - 1;
        }
    }
    return reticolo;
}                      
          


// ============================================================================
// SEZIONE 3: TOPOLOGIA E ALGORITMO A SCACCHIERA
// ============================================================================


vicini* definisci_vicini(parametri *p, int L_local) {
    int righe_totali = L_local + 2;
    vicini *lista = malloc(righe_totali * p->L * sizeof(vicini));

    // Calcolo dei vicini esclusivamente per le righe reali
    for (int j = 1; j <= L_local; j++) {
        for (int i = 0; i < p->L; i++) {

            int k = calcola_indice_lessicografico(p->L, i, j);

            // Condizioni Periodiche al Contorno (PBC) per l'asse orizzontale (X)
            int i_dx = (i + 1) % p->L; 
            int i_sx = (i - 1 + p->L) % p->L; 

            // Confini lineari per l'asse verticale (Y)
            // La periodicità globale Y è gestita dalla topologia del ring MPI
            int j_up = j + 1; 
            int j_dw = j - 1; 

            // Salvataggio indici precalcolati
            lista[k].dx = calcola_indice_lessicografico(p->L, i_dx, j);
            lista[k].sx = calcola_indice_lessicografico(p->L, i_sx, j);
            lista[k].up = calcola_indice_lessicografico(p->L, i, j_up);
            lista[k].dw = calcola_indice_lessicografico(p->L, i, j_dw);
        }
    }

    return lista;
}


void costruisci_liste_scacchiera(int L, int L_local, int rank, int **indici_rossi, int **indici_neri) {
    // In ogni dominio il 50% degli spin è rosso e il 50% è nero
    int num_per_colore = (L * L_local) / 2;

    // Allocazione memoria per i due subset indipendenti
    *indici_rossi = malloc(num_per_colore * sizeof(int));
    *indici_neri = malloc(num_per_colore * sizeof(int));

    int count_r = 0;
    int count_n = 0;

    // Analisi del solo reticolo reale locale
    for (int j = 1; j <= L_local; j++) {
        // Calcolo della coordinata Y globale assoluta per mantenere la 
        // coerenza della scacchiera attraverso i bordi dei vari domini MPI
        int y_globale = (j - 1) + (rank * L_local);

        for (int i = 0; i < L; i++) {
            int k = calcola_indice_lessicografico(L, i, j);

            // Criterio di partizionamento Red-Black basato sulle coordinate assolute
	    
            /*
	     sintassi (*puntatore)[indice]:
                 * In C, l'operatore array [] ha la precedenza sull'operatore *.
                 * Senza le tonde, il compilatore valuterebbe *(indici_rossi[count_r]), 
                 * spostando il doppio puntatore e causando un errore.
                 * Le parentesi () forzano prima la dereferenziazione per recuperare 
                 * il vero puntatore all'array nel main, e solo dopo applicano l'indice.
            */

            if ((i + y_globale) % 2 == 0) {
                (*indici_rossi)[count_r++] = k;
            } else {
                (*indici_neri)[count_n++] = k;
            }
        }
    }
}
