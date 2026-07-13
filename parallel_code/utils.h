#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

// ============================================================================
// SEZIONE 1: STRUTTURE DATI
// ============================================================================

/*
   Struttura: parametri
   --------------------
   Contiene tutti i parametri fisici e operativi necessari per la simulazione.
   Viene popolata durante la fase di lettura del file di input.
 */

typedef struct {
    // Parametri fisici
    double k_B;            // Costante di Boltzmann
    double h;              // Campo magnetico esterno
    double J;              // Costante di accoppiamento
    
    // Parametri operativi e statistiche
    double *beta_values;   // Array dei valori di beta da testare
    int num_betas;         // Dimensione dell'array beta_values
    int frame_freq;        // Intervallo di sweep per il salvataggio delle configurazioni spaziali
    int meas_freq;         // Intervallo di sweep per il campionamento delle misure termodinamiche
    int meas_sweeps;       // Numero totale di sweep dedicati alla fase di misura
    int eq_sweeps;         // Numero totale di sweep dedicati alla fase di termalizzazione
    
    // Parametri geometrici
    int L;                 // Taglia lineare corrente del reticolo
    int *L_values;         // Array dei valori di L da testare
    int num_L;             // Numero totale di taglie da testare
    int L_pow_max;         // Esponente massimo per la generazione delle taglie (es. 2^L_pow_max)
    int L_pow_min;         // Esponente minimo per la generazione delle taglie
} parametri;



// ============================================================================
// SEZIONE 2: INIZIALIZZAZIONE E LETTURA DATI
// ============================================================================

/*
   Funzione: leggi_parametri
   -------------------------
   Analizza il file di testo specificato e popola la struttura dei parametri.
   Calcola dinamicamente le taglie del reticolo utilizzando operazioni bitwise.
   Input:
   - filename: stringa contenente il percorso del file di input.
   - p: puntatore alla struttura parametri da allocare e popolare.
   Output:
   - Nessuno. Il programma termina con codice di errore se il file è inaccessibile.
 */

void leggi_parametri(const char *filename, parametri *p);

/*
   Funzione: init_rng
   ------------------
   Inizializza il generatore di numeri pseudocasuali.
   Input:
   - rank: identificativo del processo MPI corrente, utilizzato come seed
   per garantire sequenze scorrelate tra i vari processi.
   Output:
   - Nessuno.
 */
void init_rng(int rank);

// ============================================================================
// SEZIONE 3: GESTIONE INPUT/OUTPUT E SALVATAGGIO
// ============================================================================

/*
   Funzione: inizializza_file_csv
   ------------------------------
   Crea un nuovo file CSV per il salvataggio delle medie termodinamiche
   e vi inserisce l'intestazione delle colonne.
   Input:
   - filename: percorso del file CSV di destinazione.
   Output:
   - Nessuno. Il programma termina se la creazione del file fallisce.
 */
void inizializza_file_csv(const char *filename);

/*
   Funzione: salva_misura_csv
   --------------------------
   Appende una singola riga di risultati statistici al file CSV.
   Input:
   - filename: percorso del file CSV.
   - L: taglia del reticolo.
   - beta: inverso della temperatura corrente.
   - mag_media: valore atteso della magnetizzazione.
   - mag2_media: valore atteso del quadrato della magnetizzazione.
   - errore_mag: deviazione standard della magnetizzazione.
   Output:
   - Nessuno.
 */
void salva_misura_csv(const char *filename, int L, double beta, double mag_media, double mag2_media, double errore_mag);

/*
   Funzione: salva_configurazione
   ------------------------------
   Scrive l'intero stato micro-canonico del reticolo su un file di testo (formato matrice).
   Input:
   - file: puntatore al file (deve essere già aperto in modalità append).
   - reticolo: array unidimensionale contenente la configurazione corrente globale.
   - L: dimensione lineare del reticolo.
   - beta: temperatura corrente.
   - sweep: iterazione corrente della fase di misura.
   Output:
   - Nessuno.
 */
void salva_configurazione(FILE *file, int *reticolo, int L, double beta, int sweep);

// ============================================================================
// SEZIONE 4: COMUNICAZIONE MPI
// ============================================================================

/*
   Funzione: raccordo_dominio
   --------------------------
   Gestisce lo scambio delle celle fantasma (ghost cells) tra processi MPI adiacenti
   utilizzando comunicazioni bloccanti simultanee per risolvere i bordi del dominio.
   Input:
   - reticolo: array del frammento di reticolo locale (inclusivo delle righe fantasma).
   - L: dimensione lineare globale del reticolo.
   - L_local: numero di righe reali assegnate al processo corrente.
   - rank: identificativo del processo MPI corrente.
   - size: numero totale di processi allocati.
   Output:
   - Nessuno. I buffer fantasma top e bottom del reticolo vengono sovrascritti.
 */
void raccordo_dominio(int *reticolo, int L, int L_local, int rank, int size);

#endif
