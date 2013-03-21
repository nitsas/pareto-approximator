# Makefile for the experiment vs PGL's NAMOA*.
# 
# Author: Christos Nitsas
# Date:   2013
#


CC=g++
#WARNINGS=-Wall -Wextra 
WARNINGS=-Wall 
CPPLIBS=-larmadillo
INCLUDEDIR=$(PGLROOT)


all: compile

compile: Point.o main.cpp AStarDijkstra.h MultiobjectiveSpOnPmgProblem.h Hyperplane.h experiments_vs_namoa_star_common.h experiments_vs_namoa_star_utility.h
	g++ main.cpp Point.o -o vns_experiment.out -O2 $(WARNINGS) -I$(INCLUDEDIR) $(CPPLIBS)
#	g++ main.cpp Point.o -o vns_experiment.out -O2 $(WARNINGS) -I$(INCLUDEDIR) $(CPPLIBS)

debug: Point.o main.cpp AStarDijkstra.h MultiobjectiveSpOnPmgProblem.h Hyperplane.h experiments_vs_namoa_star_common.h experiments_vs_namoa_star_utility.h
	g++ main.cpp Point.o -o vns_experiment.out -O0 -g $(WARNINGS) -I$(INCLUDEDIR) $(CPPLIBS)

Point.o: ../../Point.h ../../Point.cpp ../../DifferentDimensionsException.h ../../NegativeApproximationRatioException.h ../../NotPositivePointException.h ../../NotStrictlyPositivePointException.h
	$(CC) $(WARNINGS) -c ../../Point.cpp -o $@

clean:
	rm Point.o vns_experiment.out
