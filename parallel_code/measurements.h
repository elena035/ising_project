#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include "utils.h"
#include "lattice.h"

// ============================================================================
// SEZIONE 1: OSSERVABILI TERMODINAMICHE LOCALI
// ============================================================================

/*
 Funzione: calcola_magnetizzazione_locale
 ----------------------------------------
 Calcola la somma algebrica degli spin per il frammento di reticolo
 assegnato al processo MPI corrente.
 Input:
 - reticolo: array rappresentante lo stato corrente del sistema locale.
 - L: dimensione lineare globale del reticolo.
 - L_local: numero di righe reali assegnate al processo locale.
 Output:
 - Valore di tipo double rappresentante la magnetizzazione totale del dominio locale.
 */
double calcola_magnetizzazione_locale(int *reticolo, int L, int L_local);



/*
 Funzione: calcola_energia_locale
 --------------------------------
 Calcola l'energia di interazione locale per il frammento di reticolo
 assegnato al processo MPI corrente, adottando un approccio direzionale
 per evitare il doppio conteggio dei legami.
 Input:
 - reticolo: array rappresentante lo stato corrente del sistema locale.
 - lista_vicini: mappa statica pre-calcolata dei primi vicini.
 - L: dimensione lineare globale del reticolo.
 - L_local: numero di righe reali assegnate al processo locale.
 - J: costante macroscopica di accoppiamento ferromagnetico.
 Output:
 - Valore di tipo double rappresentante l'energia termodinamica del dominio locale.
 */
double calcola_energia_locale(int *reticolo, vicini *lista_vicini, int L, int L_local, double J);

#endif
