
CC = g++

OPCIONS = -g -std=c++11 -O3 -Wall -Iinclude

all: ABC.o algorismes_comuns.o Timer.o
	$(CC) -o ABC ABC.o algorismes_comuns.o Timer.o $(OPCIONS) 

ABC: ABC.o algorismes_comuns.o Timer.o
	$(CC) -o ABC ABC.o algorismes_comuns.o Timer.o  $(OPCIONS)

ABC.o: ABC.cpp
	$(CC) -c ABC.cpp $(OPCIONS)

algorismes_comuns.o: algorismes_comuns.cpp
	$(CC) -c algorismes_comuns.cpp $(OPCIONS)

Timer.o: Timer.cc
	$(CC) -c Timer.cc $(OPCIONS)

clean:
	rm *.o ABC 
