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

echo "==================================================="
echo " Inizio Produzione: $NUM_RUNS iterazioni "
echo " Modalità: COOLING ATTIVATO (Pausa termica tra i run)"
echo "==================================================="

# 3. Ciclo FOR della produzione
for i in $(seq 1 $NUM_RUNS); do
    echo "--> Esecuzione Iterazione $i di $NUM_RUNS..."
    echo "simulazione con 1 processo"
    mpirun -np 1 ./ising $INPUT_FILE risultati_1p_def.csv
    sleep 60
    echo "simulazione con 2 processi"
    mpirun -np 2 ./ising $INPUT_FILE risultati_2p_def.csv
    sleep 60
    echo "simulazione con 4 processi"
    mpirun -np 4 ./ising $INPUT_FILE risultati_4p_def.csv
    sleep 60
    echo "simulazione con 8 processi"
    mpirun -np 8 --oversubscribe ./ising $INPUT_FILE risultati_8p_def.csv
done

# 4. Archiviazione
echo "--> Archiviazione dei dati..."
mkdir -p $DIR_OUTPUT
mv tempi_speedup_*.csv $DIR_OUTPUT/
mv risultati*.csv $DIR_OUTPUT/
mv simulazioni/ $DIR_OUTPUT/ 2>/dev/null

echo " Produzione terminata. Dati in: ./$DIR_OUTPUT/"

# 5. Avviso acustico di fine lavoro
# Se stai usando WSL su Windows, questo farà suonare la notifica di sistema
if command -v powershell.exe &> /dev/null; then
    powershell.exe -c "(New-Object System.Media.SoundPlayer 'C:\Windows\Media\tada.wav').PlaySync()"
else
    # Fallback per terminali Linux nativi o Mac (emette il 'beep' standard)
    echo -e "\a"
fi
