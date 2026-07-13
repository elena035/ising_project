#include <math.h>
#include "measurements.h"


// ============================================================================
// SEZIONE 1: OSSERVABILI TERMODINAMICHE LOCALI
// ============================================================================


double calcola_magnetizzazione_locale(int *reticolo, int L, int L_local) {
    long long somma_spin = 0;
    
    /* 
    L'iterazione è limitata strettamente alle righe reali del dominio locale
    (da j = 1 a j = L_local), escludendo le ghost cells per evitare ridondanze
    e inconsistenze durante la successiva operazione di riduzione globale (MPI_Reduce).
    */

    for (int j = 1; j <= L_local; j++) {
        for (int i = 0; i < L; i++) {
            somma_spin += reticolo[i + j * L];
        }
    }
    return (double)somma_spin;
}

double calcola_energia_locale(int *reticolo, vicini *lista_vicini, int L, int L_local, double J) {
    long long somma_energia = 0;
    
    // L'iterazione esamina unicamente i siti appartenenti al frammento locale
    for (int j = 1; j <= L_local; j++) {
        for (int i = 0; i < L; i++) {
            int k = i + j * L;
            /*
	    Al fine di prevenire il doppio conteggio delle energie di legame
            (ogni pair interaction deve essere sommata una sola volta nell'intero
            reticolo globale), si considerano unicamente i vicini in avanti
            (asse x positivo, dx) e in basso (asse y negativo, dw) rispetto al sito corrente.
	    */ 
            int interazione = reticolo[lista_vicini[k].dw] + reticolo[lista_vicini[k].dx];
            somma_energia += reticolo[k] * interazione;
        }
    }
    return (double)(-J * somma_energia);
}
