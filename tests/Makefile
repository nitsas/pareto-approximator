# Makefile for all unit tests.
# 
# Available unit tests are currently:
# - PointTest.cpp
# - Line2DTest.cpp
# - BaseProblemTest.cpp  &  NonOptimalStartingPointsProblem.h  &  
#		NonOptimalStartingPointsProblem.cpp  &  SmallGraphProblem.h  &  
#		SmallGraphProblem.cpp
# 
# Author:  Christos Nitsas
# Date:    2012
#


CC=g++
CPPFLAGS=-Wall -g
CPPLIBS=-lgtest


# Make all unit tests
all: PointTest.out Line2DTest.out BaseProblemTest.out

# Make PointTest.out
PointTest.out: PointTest.cpp ../Point.h Point.o ../DifferentDimensionsException.h ../NegativeApproximationRatioException.h ../NotPositivePointException.h
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o PointTest.cpp -o PointTest.out

# Make Line2DTest.out
Line2DTest.out: Line2DTest.cpp ../Point.h Point.o ../Line2D.h Line2D.o ../Not2DPointsException.h ../SamePointsException.h ../VerticalLineException.h ../ParallelLinesException.h
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o Line2D.o Line2DTest.cpp -o Line2DTest.out

# Make BaseProblemTest.out
BaseProblemTest.out: BaseProblemTest.o SmallGraphProblem.o NonOptimalStartingPointsProblem.o Point.o Line2D.o
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o Line2D.o NonOptimalStartingPointsProblem.o SmallGraphProblem.o BaseProblemTest.o -o BaseProblemTest.out

# Make BaseProblemTest.o
BaseProblemTest.o: BaseProblemTest.cpp ../Point.h ../PointAndSolution.h NonOptimalStartingPointsProblem.h SmallGraphProblem.h ../BaseProblem.h ../BaseProblem.cpp
	$(CC) $(CPPFLAGS) -c BaseProblemTest.cpp -o BaseProblemTest.o

# Make SmallGraphProblem.o
SmallGraphProblem.o: SmallGraphProblem.cpp SmallGraphProblem.h ../PointAndSolution.h ../BaseProblem.h ../BaseProblem.cpp ../Point.h
	$(CC) $(CPPFLAGS) -c SmallGraphProblem.cpp -o SmallGraphProblem.o

# Make NonOptimalStartingPointsProblem.o
NonOptimalStartingPointsProblem.o: NonOptimalStartingPointsProblem.cpp NonOptimalStartingPointsProblem.h ../PointAndSolution.h ../BaseProblem.h ../BaseProblem.cpp ../Point.h
	$(CC) $(CPPFLAGS) -c NonOptimalStartingPointsProblem.cpp -o NonOptimalStartingPointsProblem.o


# Make Point.o
Point.o: ../Point.h ../Point.cpp ../DifferentDimensionsException.h ../NegativeApproximationRatioException.h ../NotPositivePointException.h
	$(CC) $(CPPFLAGS) -c ../Point.cpp -o Point.o

# Make Line2D.o
Line2D.o: ../Line2D.h ../Line2D.cpp ../Point.h Point.o ../Not2DPointsException.h ../SamePointsException.h ../VerticalLineException.h ../ParallelLinesException.h
	$(CC) $(CPPFLAGS) -c ../Line2D.cpp -o Line2D.o

# Clean object files and executables
clean: 
	rm -f Point.o Line2D.o BaseProblemTest.o SmallGraphProblem.o NonOptimalStartingPointsProblem.o PointTest.out Line2DTest.out BaseProblemTest.out

