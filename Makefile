CC = g++
CFLAGS= -g -Wall -std=c++0x

default: compile

compile: Graph.o Aco.o Main.o
	$(CC) $? -o aco

Main.o: Main.cpp
	$(CC) $(CFLAGS) -c $< -o $@

Graph.o: Graph.cpp Graph.h
	$(CC) $(CFLAGS) -c $< -o $@

Aco.o: ACO.cpp ACO.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) *.o *~ aco
