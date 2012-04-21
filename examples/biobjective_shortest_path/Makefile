# Makefile for the biobjective_shortest_path example.
#
# Author:  Christos Nitsas
# Date:    2012
# 


CC=g++
CPPFLAGS=-Wall -Wextra -Werror
DEBUGFLAGS=-g


# Link everything and make example1
bosp_example.out: main.o ../../Point.h Point.o ../../Line2D.h Line2D.o ../../BaseProblem.h ../../BaseProblem.cpp RandomGraphProblem.h RandomGraphProblem.o ../../PointAndSolution.h
	$(CC) $(CPPFLAGS) $(DEBUGFLAGS) Point.o Line2D.o RandomGraphProblem.o main.o -o $@

# Compile main.cpp
main.o: main.cpp ../../Point.h ../../Line2D.h ../../BaseProblem.h ../../BaseProblem.cpp RandomGraphProblem.h ../../PointAndSolution.h
	$(CC) $(CPPFLAGS) $(DEBUGFLAGS) -c main.cpp -o main.o

# Make RandomGraphProblem.o
RandomGraphProblem.o: RandomGraphProblem.h RandomGraphProblem.cpp ../../PointAndSolution.h ../../BaseProblem.h ../../BaseProblem.cpp ../../Point.h ../../Line2D.h
	$(CC) $(CPPFLAGS) $(DEBUGFLAGS) -c RandomGraphProblem.cpp -o RandomGraphProblem.o


# Make Point.o
Point.o: ../../Point.h ../../Point.cpp ../../DifferentDimensionsException.h ../../NonExistentCoordinateException.h ../../NegativeApproximationRatioException.h ../../NotPositivePointException.h
	$(CC) $(CPPFLAGS) $(DEBUGFLAGS) -c ../../Point.cpp -o Point.o

# Make Line2D.o
Line2D.o: ../../Line2D.h ../../Line2D.cpp ../../Point.h ../../Not2DPointsException.h ../../SamePointsException.h ../../VerticalLineException.h ../../ParallelLinesException.h
	$(CC) $(CPPFLAGS) $(DEBUGFLAGS) -c ../../Line2D.cpp -o Line2D.o

# Make debug version


# Clean object files and executables
clean: 
	rm -f main.o RandomGraphProblem.o Point.o Line2D.o bosp_example.out

