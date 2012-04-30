# Makefile for the biobjective_shortest_path example.
#
# Author:  Christos Nitsas
# Date:    2012
# 


CC=g++
CPPFLAGS=-Wall -Wextra -Werror
DEBUGFLAGS=-g


# Link everything and make example1
bosp_example.out: main.o ../../Point.h Point.o ../../Hyperplane.h Hyperplane.o ../../BaseProblem.h ../../BaseProblem.cpp RandomGraphProblem.h RandomGraphProblem.o ../../PointAndSolution.h ../../PointAndSolution.cpp
	$(CC) $(CPPFLAGS) $(DEBUGFLAGS) Point.o Hyperplane.o RandomGraphProblem.o main.o -o $@

# Compile main.cpp
main.o: main.cpp ../../BaseProblem.h ../../BaseProblem.cpp RandomGraphProblem.h ../../PointAndSolution.h ../../PointAndSolution.cpp
	$(CC) $(CPPFLAGS) $(DEBUGFLAGS) -c main.cpp -o main.o

# Make RandomGraphProblem.o
RandomGraphProblem.o: RandomGraphProblem.h RandomGraphProblem.cpp ../../PointAndSolution.h ../../PointAndSolution.cpp ../../BaseProblem.h ../../BaseProblem.cpp ../../Point.h 
	$(CC) $(CPPFLAGS) $(DEBUGFLAGS) -c RandomGraphProblem.cpp -o RandomGraphProblem.o


# Make Point.o
Point.o: ../../Point.h ../../Point.cpp ../../DifferentDimensionsException.h ../../NonExistentCoordinateException.h ../../NegativeApproximationRatioException.h ../../NotPositivePointException.h
	$(CC) $(CPPFLAGS) $(DEBUGFLAGS) -c ../../Point.cpp -o Point.o

# Make Hyperplane.o
Hyperplane.o: ../../Hyperplane.h ../../Hyperplane.cpp ../../Point.h Point.o ../../DifferentDimensionsException.h ../../SamePointsException.h ../../Not2DPointsException.h ../../Not2DHyperplanesException.h ../../NonExistentCoefficientException.h ../../ParallelHyperplanesException.h
	$(CC) $(CPPFLAGS) -c ../../Hyperplane.cpp -o Hyperplane.o


# Clean object files and executables
clean: 
	rm -f main.o RandomGraphProblem.o Point.o Hyperplane.o bosp_example.out

