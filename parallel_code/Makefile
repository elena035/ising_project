# Compilatore e flag
CC = mpicc
CFLAGS = -Wall -O3 -ffast-math -march=native
LDFLAGS = -lm

# Nome dell'eseguibile finale
TARGET = ising

# Trova tutti i file .c nella cartella
SRCS = main.c utils.c lattice.c metropolis.c measurements.c validazione.c
OBJS = $(SRCS:.c=.o)

# Regola principale
all: $(TARGET)

# Come costruire l'eseguibile unendo gli oggetti
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Come compilare i singoli file .c in file .o
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Regola per pulire la cartella dai file compilati
clean:
	rm -f $(OBJS) $(TARGET)
