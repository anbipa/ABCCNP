
CC = g++

OPCIONS = -g -std=c++11 -O3

all: MISGreedy.o bestSwapApproach.o ComponentBased.o ComponentBasedRanked.o VCoverOrdered.o VCoverShuffled.o NonDetVertexCover.o Greedy3d.o CompBasedGreedy3d.o prova.o algorismes_comuns.o
	$(CC) -o MISGreedy MISGreedy.o algorismes_comuns.o $(OPCIONS) 
	$(CC) -o bestSwapApproach bestSwapApproach.o algorismes_comuns.o $(OPCIONS)
	$(CC) -o ComponentBased ComponentBased.o algorismes_comuns.o $(OPCIONS)
	$(CC) -o ComponentBasedRanked ComponentBasedRanked.o algorismes_comuns.o $(OPCIONS)
	$(CC) -o VCoverOrdered VCoverOrdered.o algorismes_comuns.o $(OPCIONS)
	$(CC) -o VCoverShuffled VCoverShuffled.o algorismes_comuns.o $(OPCIONS)
	$(CC) -o NonDetVertexCover NonDetVertexCover.o algorismes_comuns.o $(OPCIONS)
	$(CC) -o Greedy3d Greedy3d.o algorismes_comuns.o $(OPCIONS)
	$(CC) -o CompBasedGreedy3d CompBasedGreedy3d.o algorismes_comuns.o $(OPCIONS)
	$(CC) -o prova prova.o algorismes_comuns.o $(OPCIONS)


MISGreedy: MISGreedy.o algorismes_comuns.o
	$(CC) -o MISGreedy MISGreedy.o algorismes_comuns.o  $(OPCIONS)

bestSwapApproach: bestSwapApproach.o algorismes_comuns.o
	$(CC) -o bestSwapApproach bestSwapApproach.o algorismes_comuns.o  $(OPCIONS)

ComponentBased: ComponentBased.o algorismes_comuns.o
	$(CC) -o ComponentBased ComponentBased.o algorismes_comuns.o  $(OPCIONS)

ComponentBasedRanked: ComponentBasedRanked.o algorismes_comuns.o
	$(CC) -o ComponentBasedRanked ComponentBasedRanked.o algorismes_comuns.o  $(OPCIONS)

VCoverOrdered: VCoverOrdered.o algorismes_comuns.o
	$(CC) -o VCoverOrdered VCoverOrdered.o algorismes_comuns.o  $(OPCIONS)

VCoverShuffled: VCoverShuffled.o algorismes_comuns.o
	$(CC) -o VCoverShuffled VCoverShuffled.o algorismes_comuns.o  $(OPCIONS)

NonDetVertexCover: NonDetVertexCover.o algorismes_comuns.o
	$(CC) -o NonDetVertexCover NonDetVertexCover.o algorismes_comuns.o  $(OPCIONS)

Greedy3d: Greedy3d.o algorismes_comuns.o
	$(CC) -o Greedy3d Greedy3d.o algorismes_comuns.o  $(OPCIONS)

CompBasedGreedy3d: CompBasedGreedy3d.o algorismes_comuns.o
	$(CC) -o CompBasedGreedy3d CompBasedGreedy3d.o algorismes_comuns.o  $(OPCIONS)

prova: prova.o algorismes_comuns.o
	$(CC) -o prova prova.o algorismes_comuns.o  $(OPCIONS)



MISGreedy.o: MISGreedy.cpp
	$(CC) -c MISGreedy.cpp $(OPCIONS)

bestSwapApproach.o: bestSwapApproach.cpp
	$(CC) -c bestSwapApproach.cpp $(OPCIONS)

ComponentBased.o: ComponentBased.cpp
	$(CC) -c ComponentBased.cpp $(OPCIONS)

ComponentBasedRanked.o: ComponentBasedRanked.cpp
	$(CC) -c ComponentBasedRanked.cpp $(OPCIONS)

VCoverOrdered.o: VCoverOrdered.cpp
	$(CC) -c VCoverOrdered.cpp $(OPCIONS)

VCoverShuffled.o: VCoverShuffled.cpp
	$(CC) -c VCoverShuffled.cpp $(OPCIONS)

NonDetVertexCover.o: NonDetVertexCover.cpp
	$(CC) -c NonDetVertexCover.cpp $(OPCIONS)

Greedy3d.o: Greedy3d.cpp
	$(CC) -c Greedy3d.cpp $(OPCIONS)

CompBasedGreedy3d.o: CompBasedGreedy3d.cpp
	$(CC) -c CompBasedGreedy3d.cpp $(OPCIONS)

algorismes_comuns.o: algorismes_comuns.cpp
	$(CC) -c algorismes_comuns.cpp $(OPCIONS)

prova.o: prova.cpp
	$(CC) -c prova.cpp $(OPCIONS)



clean:
	rm *.o MISGreedy bestSwapApproach ComponentBased ComponentBasedRanked VCoverOrdered VCoverShuffled NonDetVertexCover Greedy3d CompBasedGreedy3d prova
