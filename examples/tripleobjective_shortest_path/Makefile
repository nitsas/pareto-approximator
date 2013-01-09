# Makefile for the tripleobjective_shortest_path example.
#
# Author:  Christos Nitsas
# Date:    2012
# 


CC=g++
CPPFLAGS=-Wall -Wextra -Werror -g -O2
CPPLIBS=-larmadillo


# Link everything and make tosp_example.out
tosp_example.out: main.o ../../Point.h Point.o ../../BaseProblem.h ../../BaseProblem.cpp RandomGraphProblem.h RandomGraphProblem.o ../../PointAndSolution.h ../../PointAndSolution.cpp FloodVisitor.o FloodVisitor.h
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o RandomGraphProblem.o main.o FloodVisitor.o -o $@

# Compile main.cpp
main.o: main.cpp ../../utility.h ../../utility.cpp ../../BaseProblem.h ../../BaseProblem.cpp ../../Facet.h ../../Facet.cpp RandomGraphProblem.h ../../PointAndSolution.h ../../PointAndSolution.cpp tripleobjective_shortest_path_example_common.h
	$(CC) $(CPPFLAGS) -c main.cpp -o $@

# Make RandomGraphProblem.o
RandomGraphProblem.o: RandomGraphProblem.h RandomGraphProblem.cpp ../../PointAndSolution.h ../../PointAndSolution.cpp ../../Facet.h ../../Facet.cpp ../../BaseProblem.h ../../BaseProblem.cpp ../../utility.h ../../utility.cpp ../../Point.h tripleobjective_shortest_path_example_common.h
	$(CC) $(CPPFLAGS) -c RandomGraphProblem.cpp -o $@

# Make FloodVisitor.o
FloodVisitor.o: FloodVisitor.h FloodVisitor.cpp tripleobjective_shortest_path_example_common.h ../../Point.h ../../NonDominatedSet.h ../../NonDominatedSet.cpp
	$(CC) $(CPPFLAGS) -c FloodVisitor.cpp -o $@


# Make Point.o
Point.o: ../../Point.h ../../Point.cpp ../../NullObjectException.h ../../DifferentDimensionsException.h ../../NonExistentCoordinateException.h ../../NegativeApproximationRatioException.h ../../NotPositivePointException.h ../../NotStrictlyPositivePointException.h
	$(CC) $(CPPFLAGS) -c ../../Point.cpp -o $@


# Clean object files and executables
clean: 
	rm -f main.o RandomGraphProblem.o FloodVisitor.o Point.o tosp_example.out

