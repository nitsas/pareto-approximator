# Makefile for the biobjective_shortest_path example.
#
# Author:  Christos Nitsas
# Date:    2012
# 


CC=g++
CPPFLAGS=-Wall -Wextra -Werror -g
CPPLIBS=-larmadillo


# Link everything and make bosp_example.out
bosp_example.out: main.o ../../Point.h Point.o ../../Hyperplane.h Hyperplane.o ../../BaseProblem.h ../../BaseProblem.cpp RandomGraphProblem.h RandomGraphProblem.o ../../PointAndSolution.h ../../PointAndSolution.cpp FloodVisitor.o FloodVisitor.h
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o Hyperplane.o RandomGraphProblem.o main.o FloodVisitor.o -o $@

# Compile main.cpp
main.o: main.cpp ../../BaseProblem.h ../../BaseProblem.cpp RandomGraphProblem.h ../../PointAndSolution.h ../../PointAndSolution.cpp biobjective_shortest_path_example_common.h
	$(CC) $(CPPFLAGS) -c main.cpp -o $@

# Make RandomGraphProblem.o
RandomGraphProblem.o: RandomGraphProblem.h RandomGraphProblem.cpp ../../PointAndSolution.h ../../PointAndSolution.cpp ../../BaseProblem.h ../../BaseProblem.cpp ../../Point.h biobjective_shortest_path_example_common.h
	$(CC) $(CPPFLAGS) -c RandomGraphProblem.cpp -o $@

# Make FloodVisitor.o
FloodVisitor.o: FloodVisitor.h FloodVisitor.cpp biobjective_shortest_path_example_common.h ../../Point.h ../../NonDominatedSet.h ../../NonDominatedSet.cpp
	$(CC) $(CPPFLAGS) -c FloodVisitor.cpp -o $@


# Make Point.o
Point.o: ../../Point.h ../../Point.cpp ../../DifferentDimensionsException.h ../../NonExistentCoordinateException.h ../../NegativeApproximationRatioException.h ../../NotPositivePointException.h
	$(CC) $(CPPFLAGS) -c ../../Point.cpp -o $@

# Make Hyperplane.o
Hyperplane.o: ../../Hyperplane.h ../../Hyperplane.cpp ../../Point.h Point.o ../../DifferentDimensionsException.h ../../NonExistentCoefficientException.h ../../PointIsOnTheHyperplaneException.h
	$(CC) $(CPPFLAGS) -c ../../Hyperplane.cpp -o $@


# Clean object files and executables
clean: 
	rm -f main.o RandomGraphProblem.o FloodVisitor.o Point.o Hyperplane.o bosp_example.out

