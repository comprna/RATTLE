EXEC=rattle
CC=g++
CFLAGS=-Wall -Wextra -std=c++14 -O3 -pthread

all: $(EXEC)

$(EXEC): main.cpp fasta.o cluster.o utils.o kmer.o similarity.o correct.o
	$(CC) -o $(EXEC) $(CFLAGS) main.cpp fasta.o cluster.o utils.o kmer.o similarity.o correct.o

utils.o: utils.hpp utils.cpp
	$(CC) -c $(CFLAGS) utils.cpp

fasta.o: fasta.hpp fasta.cpp
	$(CC) -c $(CFLAGS) fasta.cpp

kmer.o: kmer.hpp kmer.cpp
	$(CC) -c $(CFLAGS) kmer.cpp

similarity.o: similarity.hpp similarity.cpp
	$(CC) -c $(CFLAGS) similarity.cpp

cluster.o: cluster.hpp cluster.cpp kmer.o similarity.o utils.o fasta.o
	$(CC) -c $(CFLAGS) cluster.cpp kmer.cpp similarity.cpp utils.cpp fasta.cpp

correct.o: correct.hpp correct.cpp utils.o fasta.o
	$(CC) -c $(CFLAGS) correct.cpp utils.cpp fasta.cpp

clean: 
	rm -f *.o
	rm -f $(EXEC)
	rm -f *log.txt