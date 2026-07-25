#!/bin/bash

# --- CONFIGURAZIONE ---
NUM_RUNS=1
INPUT_FILE="input.txt" 
# ----------------------

# 1. Creiamo il nome della cartella basato sulla data e ora attuale
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
DIR_OUTPUT="produzione_$TIMESTAMP"

echo "=== Compilazione ==="
make clean
make

mkdir -p vecchi_test
mv tempi_speedup_*.csv vecchi_test/ 2>/dev/null
mv simulazioni/ vecchi_test/ 2>/dev/null
mv risultati_*.csv vecchi_test/ 2>/dev/null


mkdir -p "$DIR_OUTPUT"


echo "==================================================="
echo " Inizio Produzione: $NUM_RUNS iterazioni "
echo " Modalità: COOLING ATTIVATO (Pausa termica tra i run)"
echo "==================================================="

# 3. Ciclo FOR della produzione
for i in $(seq 1 $NUM_RUNS); do
    echo "--> Esecuzione Iterazione $i di $NUM_RUNS..."
    
    echo "simulazione con 1 processo"
    mpirun -np 1 ./ising $INPUT_FILE risultati_1p.csv
    sleep 60
    
    echo "simulazione con 2 processi"
    mpirun -np 2 ./ising $INPUT_FILE risultati_2p.csv
    sleep 60
    
    echo "simulazione con 3 processi"
    mpirun -np 3 ./ising $INPUT_FILE risultati_3p.csv
    sleep 60
    
    echo "simulazione con 4 processi"
    mpirun -np 4 ./ising $INPUT_FILE risultati_4p.csv
    sleep 60
    
    echo "simulazione con 5 processi"
    mpirun -np 5 ./ising $INPUT_FILE risultati_5p.csv

    echo "simulazione con 6 processi"
    mpirun -np 6 ./ising $INPUT_FILE risultati_6p.csv   
done

echo "==================================================="
echo "Spostamento dei risultati in $DIR_OUTPUT..."
# SPOSTAMENTO DEI FILE GENERATI NELLA CARTELLA DI PRODUZIONE
mv risultati_*.csv "$DIR_OUTPUT"/ 2>/dev/null
mv tempi_speedup_*.csv "$DIR_OUTPUT"/ 2>/dev/null
mv simulazioni/ "$DIR_OUTPUT"/ 2>/dev/null

echo "Produzione terminata con successo!"
