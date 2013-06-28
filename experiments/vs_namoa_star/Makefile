# Makefile for the experiment vs PGL's NAMOA*.
# 
# Author: Christos Nitsas
# Date:   2013
#


CC=g++
# the ${HOME}/usr/lib64 and ${HOME}/usr/include folders below are for the armadillo library
CPPFLAGS=-Wall -L${HOME}/usr/lib64 -I${HOME}/usr/include -Wno-deprecated -Wno-unused-but-set-variable
CPPLIBS=-larmadillo
INCLUDEDIR=$(PGLROOT)


all: compile

compile: Point.o main.cpp AStarDijkstra.h MultiobjectiveSpOnPmgProblem.h Hyperplane.h experiments_vs_namoa_star_common.h experiments_vs_namoa_star_utility.h
	g++ main.cpp Point.o -o vns_experiment.out -O2 $(CPPFLAGS) -I$(INCLUDEDIR) $(CPPLIBS)

debug: Point.o main.cpp AStarDijkstra.h MultiobjectiveSpOnPmgProblem.h Hyperplane.h experiments_vs_namoa_star_common.h experiments_vs_namoa_star_utility.h
	g++ main.cpp Point.o -o vns_experiment.out -O0 -g $(CPPFLAGS) -I$(INCLUDEDIR) $(CPPLIBS)

Point.o: ../../Point.h ../../Point.cpp ../../DifferentDimensionsException.h ../../NegativeApproximationRatioException.h ../../NotPositivePointException.h ../../NotStrictlyPositivePointException.h
	$(CC) $(CPPFLAGS) -c ../../Point.cpp -o $@

clean:
	rm Point.o vns_experiment.out
