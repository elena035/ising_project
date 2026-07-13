// lattice.h
#ifndef LATTICE_H
#define LATTICE_H

#include "utils.h"

// ============================================================================
// SEZIONE 1: STRUTTURE DATI
// ============================================================================


/*
 Struttura: coordinate
 ---------------------
 Rappresenta la posizione spaziale di un sito nel reticolo bidimensionale.
 */
typedef struct { 
    int x; // Indice di colonna (asse orizzontale)
    int y; // Indice di riga (asse verticale)
} coordinate;

/*
 Struttura: vicini
 -----------------
 Contiene gli indici lessicografici dei quattro siti primi vicini 
 di un dato nodo nel reticolo.
 */
typedef struct { 
    int up; // Indice del vicino superiore
    int dw; // Indice del vicino inferiore
    int dx; // Indice del vicino destro
    int sx; // Indice del vicino sinistro
} vicini;



// ============================================================================
// SEZIONE 2: GEOMETRIA E MAPPING
// ============================================================================


/*
 Funzione: calcola_indice_lessicografico
 ---------------------------------------
 Converte le coordinate cartesiane (x, y) in un indice monodimensionale (array 1D).
 Input:
 - L: dimensione lineare del reticolo globale.
 - i: coordinata x (colonna).
 - j: coordinata y (riga).
 Output:
 - Intero rappresentante l'indice lessicografico (k).
 */
int calcola_indice_lessicografico(int L, int i, int j);

/*
 Funzione: calcola_coordinate_cartesiane
 ---------------------------------------
 Operazione inversa: converte un indice monodimensionale nelle rispettive
 coordinate cartesiane (x, y).
 Input:
 - k: indice lessicografico nell'array 1D.
 - L: dimensione lineare del reticolo globale.
 Output:
 - Struttura 'coordinate' contenente i valori (x, y).
 */
coordinate calcola_coordinate_cartesiane(int k, int L);



// ============================================================================
// SEZIONE 3: INIZIALIZZAZIONE E TOPOLOGIA
// ============================================================================

/*
 Funzione: costruisci_reticolo
 -----------------------------
 Alloca la memoria per il frammento di reticolo locale assegnato al processo MPI,
 includendo le righe fantasma, e lo inizializza con una configurazione random.
 Input:
 - p: puntatore alla struttura parametri.
 - L_local: numero di righe reali assegnate al processo corrente.
 Output:
 - Puntatore all'array di interi rappresentante il reticolo locale inizializzato.
 */
int* costruisci_reticolo(parametri *p, int L_local);


/*
 Funzione: definisci_vicini
 --------------------------
 Costruisce una mappa statica dei primi vicini per ogni sito del reticolo locale.
 Applica le condizioni periodiche al contorno (PBC) sull'asse X, mentre lascia
 l'asse Y lineare per delegare il raccordo a MPI.
 Input:
 - p: puntatore alla struttura parametri.
 - L_local: numero di righe reali assegnate al processo corrente.
 Output:
 - Puntatore all'array di strutture 'vicini' precalcolato.
 */
vicini* definisci_vicini(parametri *p, int L_local);


/*
 Funzione: costruisci_liste_scacchiera
 -------------------------------------
 Suddivide gli indici del dominio locale in due set indipendenti (rosso e nero)
 per consentire l'aggiornamento parallelo senza rompere il bilancio dettagliato (Red-Black updating).
 Input:
 - L: dimensione lineare globale del reticolo.
 - L_local: numero di righe reali assegnate al processo corrente.
 - rank: identificativo del processo MPI corrente.
 - indici_rossi: doppio puntatore per l'allocazione dell'array degli indici rossi.
 - indici_neri: doppio puntatore per l'allocazione dell'array degli indici neri.
 Output:
 - Nessuno (gli array vengono allocati e popolati per riferimento).
 */
void costruisci_liste_scacchiera(int L, int L_local, int rank, int **indici_rossi, int **indici_neri);

#endif
