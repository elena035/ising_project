#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "utils.h"
#include "lattice.h"

// ============================================================================
// SEZIONE 1: FISICA DEL MODELLO E PROBABILITA'
// ============================================================================


/*
Funzione: init_probabilities
----------------------------
Pre-calcola e memorizza le probabilità di accettazione di Metropolis per i soli
valori di Delta(E) positivi (ovvero 4 e 8), evitando milioni di chiamate ridondanti 
alla funzione esponenziale exp() durante il ciclo principale.
Input:
- beta: inverso della temperatura corrente.
- J: costante di accoppiamento.
- probabilities: array pre-allocato (dimensione minima 9) per salvare i valori.
Output:
- Nessuno. L'array viene popolato per riferimento.
*/
void init_probabilities(double beta, double J, double *probabilities);

/*
Funzione: calcola_delta_E
-------------------------
Calcola la variazione di energia locale conseguente alla potenziale inversione 
di un singolo spin, basandosi sull'interazione con i suoi 4 primi vicini.
Input:
- reticolo: array rappresentante lo stato corrente del sistema.
- lista_vicini: mappa statica dei primi vicini pre-calcolata.
- k: indice lessicografico dello spin da valutare.
Output:
- Intero che rappresenta la variazione di energia Delta(E).
*/
int calcola_delta_E(int *reticolo, vicini *lista_vicini, int k);



// ============================================================================
// SEZIONE 2: DINAMICA E AGGIORNAMENTO
// ============================================================================

/*
Funzione: metropolis_mezzo_sweep
--------------------------------
Applica l'algoritmo di Metropolis a un sottoinsieme di spin spazialmente
indipendenti (solo colore rosso o solo colore nero) scelti in ordine sequenziale.
Input:
- reticolo: array rappresentante lo stato corrente del sistema.
- lista_vicini: mappa statica dei primi vicini.
- probabilities: array delle probabilità esponenziali pre-calcolate.
- indici_colore: array contenente gli indici degli spin di un singolo colore.
- num_indici: dimensione dell'array indici_colore.
Output:
- Nessuno. Il reticolo viene aggiornato in-place.
*/
void metropolis_mezzo_sweep(int *reticolo, vicini *lista_vicini, double *probabilities, int *indici_colore, int num_indici);


/*
Funzione: metropolis_sweep_completo
-----------------------------------
Orchestra un intero step temporale (sweep) per il processo MPI corrente.
Utilizza il partizionamento a scacchiera per alternare cicli di calcolo puro
(rosso/nero) a cicli di comunicazione MPI necessari per aggiornare le ghost cells.
Input:
- reticolo: array rappresentante lo stato corrente del sistema locale.
- lista_vicini: mappa statica dei primi vicini.
- probabilities: array delle probabilità esponenziali pre-calcolate.
- indici_rossi: array con gli indici del sotto-reticolo rosso.
- indici_neri: array con gli indici del sotto-reticolo nero.
- L: dimensione lineare globale del reticolo.
- L_local: numero di righe reali assegnate al processo.
- rank: identificativo del processo MPI corrente.
- size: numero totale di processi MPI.
- t_comp: puntatore all'accumulatore per misurare il tempo di puro calcolo.
- t_comm: puntatore all'accumulatore per misurare il tempo di comunicazione MPI.
Output:
- Nessuno. Il reticolo locale e gli accumulatori di tempo vengono aggiornati.
*/
void metropolis_sweep_completo(int *reticolo, vicini *lista_vicini, double *probabilities, int *indici_rossi, int *indici_neri, int L, int L_local, int rank, int size, double *t_comp, double *t_comm);


#endif
