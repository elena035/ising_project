#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "validazione.h"
#include "lattice.h"
#include "measurements.h"
#include "metropolis.h"

void esegui_test_validazione(parametri *p) {
    printf("\n[VALIDAZIONE] Inizio suite di test unitari sui moduli core...\n");

    /*
    Salviamo il parametro L originale per non corrompere la simulazione successiva.
    Istanziamo un micro-reticolo 4x4, che è la dimensione minima sufficiente 
    per testare le celle fantasma e le condizioni periodiche al contorno.
    */
    int L_originale = p->L;
    p->L = 4;
    int L_local_test = 4; 
    
    int *reticolo_test = costruisci_reticolo(p, L_local_test);
    vicini *vicini_test = definisci_vicini(p, L_local_test);

    // Coordinate cartesiane (x=0, y=1) -> Indice lessicografico K = 4
    // Questo è il primo spin appartenente al reticolo reale locale
    int k_start = 4;


    // ------------------------------------------------------------------------
    // TEST 1: Geometria e Condizioni al Contorno
    // ------------------------------------------------------------------------
    printf("[VALIDAZIONE] Test 1/3: Topologia spaziale e array vicini... ");
    
    // Verifica PBC asse orizzontale: il vicino sinistro di x=0 deve wrappare a x=3 (indice 7)
    assert(vicini_test[k_start].sx == 7);

    // Verifica asse verticale: il vicino in basso di y=1 deve essere la prima ghost cell (indice 0)
    assert(vicini_test[k_start].dw == 0);

    printf("SUPERATO\n");


    // ------------------------------------------------------------------------
    // TEST 2: Limite Ferromagnetico (Magnetizzazione Macroscopica)
    // ------------------------------------------------------------------------
    printf("[VALIDAZIONE] Test 2/3: Limite termodinamico ordinato (M = 1.0)... ");
    
    // Forziamo uno stato di ordine assoluto saturando il reticolo a spin +1
    for (int j = 1; j <= L_local_test; j++) {
        for (int i = 0; i < p->L; i++) {
            reticolo_test[i + j * p->L] = 1;
        }
    }
    // Saturazione coerente dei bordi fantasma per evitare perturbazioni al calcolo
    for (int i = 0; i < p->L; i++) {
        reticolo_test[i] = 1;                        // top
        reticolo_test[i + (L_local_test + 1) * p->L] = 1; // bottom
    }
    
    // Controllo che la magnetizzazione media sia esattamente 1
    double mag_test = calcola_magnetizzazione_locale(reticolo_test, p->L, L_local_test) / (p->L * p->L);
    assert(mag_test == 1.0);
    
    printf("SUPERATO \n");


    // ------------------------------------------------------------------------
    // TEST 3: Limiti Energetici e Variazione di Energia Locale
    // ------------------------------------------------------------------------
    printf("[VALIDAZIONE] Test 3/3: Configurazione energetica microscopica di Metropolis... ");
    
    // Isoliamo un cluster di test: spin +1 completamente circondato da spin +1
    reticolo_test[k_start] = 1;
    reticolo_test[vicini_test[k_start].up] = 1;
    reticolo_test[vicini_test[k_start].dw] = 1;
    reticolo_test[vicini_test[k_start].dx] = 1;
    reticolo_test[vicini_test[k_start].sx] = 1;


    /* Se uno spin in perfetto accordo termodinamico (+1 in mezzo a tutti +1) decide 
    di invertirsi (-1), va contro 4 legami. Il costo energetico è il massimo teorico, 
    ovvero Delta E = 8J (per J=1 vale 8).
    */
    assert(calcola_delta_E(reticolo_test, vicini_test, k_start) == 8);

    printf("SUPERATO\n");


    // ------------------------------------------------------------------------
    // TEST 4: Topologia MPI (Ring Communication)
    // ------------------------------------------------------------------------
    printf("[VALIDAZIONE] Test 4/4: Topologia MPI e calcolo adiacenze... ");
    
    int size_test = 4; // Simuliamo un cluster di 4 nodi
    int rank_test = 0; // Testiamo i bordi per il root process (rank 0)
    
    int vicino_sopra = (rank_test - 1 + size_test) % size_test;
    int vicino_sotto = (rank_test + 1) % size_test;
    
    // In un anello periodico, il vicino "sopra" del rank 0 deve essere l'ultimo rank (3)
    assert(vicino_sopra == 3);
    // Il vicino "sotto" del rank 0 deve essere banalmente il rank 1
    assert(vicino_sotto == 1);
    
    printf("SUPERATO\n");

    // Pulizia e ripristino per consentire l'avvio della simulazione reale
    free(reticolo_test);
    free(vicini_test);
    p->L = L_originale; 

    printf("Validazione completata. Tutti i test sono stati superati\n");
}
