EXEC=rattle
CC=g++
CFLAGS=-Wall -Wextra -std=c++14 -O3 -pthread -I./seqan/include/

all: $(EXEC)

$(EXEC): main.cpp fasta.cpp utils.cpp kmer.cpp similarity.cpp cluster.cpp
	$(CC) -o $(EXEC) $(CFLAGS) main.cpp fasta.cpp utils.cpp kmer.cpp similarity.cpp cluster.cpp

clean: 
	rm -f *.o
	rm -f $(EXEC)
	rm -f *log.txt
