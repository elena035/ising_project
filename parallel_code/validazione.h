#ifndef VALIDAZIONE_H
#define VALIDAZIONE_H

#include "utils.h"

// ============================================================================
// SEZIONE 1: TEST DI VALIDAZIONE (UNIT TESTING)
// ============================================================================

/*
Funzione: esegui_test_validazione
---------------------------------
Esegue una suite di test unitari automatizzati per verificare la correttezza 
matematica e strutturale dei moduli core del programma: mapping spaziale, 
osservabili termodinamiche macroscopiche ed energia microscopica.

Input:
- p: puntatore alla struttura parametri (per eseguire il backup temporaneo di L).

Output:
- Nessuno. Stampa a terminale i log di progresso ed esito dei test. Interrompe 
l'esecuzione tramite assert() qualora si verifichino anomalie o bug computazionali.
*/
void esegui_test_validazione(parametri *p);

#endif
