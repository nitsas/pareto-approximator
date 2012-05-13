# Makefile for all unit tests.
# 
# Available unit tests are currently:
# - PointTest.cpp
# - Line2DTest.cpp
# - HyperplaneTest.cpp
# - BaseProblemTest.cpp  &  NonOptimalStartingPointsProblem.h  &  
#		NonOptimalStartingPointsProblem.cpp  &  SmallBiobjectiveSPProblem.h  &  
#		SmallBiobjectiveSPProblem.cpp $ SmallTripleobjectiveSPProblem.h & 
#   SmallTrimpleobjectiveSPProblem.cpp
# 
# Author:  Christos Nitsas
# Date:    2012
#


CC=g++
CPPFLAGS=-Wall -Wextra -Werror -g
CPPLIBS=-lgtest -larmadillo


# Make all unit tests
all: PointTest.out Line2DTest.out HyperplaneTest.out BaseProblemTest.out

# Make PointTest.out
PointTest.out: PointTest.cpp ../Point.h Point.o ../DifferentDimensionsException.h ../NegativeApproximationRatioException.h ../NotPositivePointException.h ../NonExistentCoordinateException.h
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o PointTest.cpp -o PointTest.out

# Make Line2DTest.out
Line2DTest.out: Line2DTest.cpp ../Point.h Point.o ../Line2D.h Line2D.o ../Not2DPointsException.h ../SamePointsException.h ../VerticalLineException.h ../ParallelLinesException.h
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o Line2D.o Line2DTest.cpp -o Line2DTest.out

# Make HyperplaneTest.out
HyperplaneTest.out: HyperplaneTest.cpp ../Point.h Point.o ../Hyperplane.h Hyperplane.o ../NonExistentCoefficientException.h
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o Hyperplane.o HyperplaneTest.cpp -o HyperplaneTest.out

# Make BaseProblemTest.out
BaseProblemTest.out: BaseProblemTest.o SmallBiobjectiveSPProblem.o SmallTripleobjectiveSPProblem.o NonOptimalStartingPointsProblem.o Point.o Hyperplane.o
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o Hyperplane.o NonOptimalStartingPointsProblem.o SmallBiobjectiveSPProblem.o SmallTripleobjectiveSPProblem.o BaseProblemTest.o -o BaseProblemTest.out

# Make BaseProblemTest.o
BaseProblemTest.o: BaseProblemTest.cpp ../Point.h ../PointAndSolution.h ../PointAndSolution.cpp NonOptimalStartingPointsProblem.h SmallBiobjectiveSPProblem.h SmallTripleobjectiveSPProblem.h ../BaseProblem.h ../BaseProblem.cpp
	$(CC) $(CPPFLAGS) -c BaseProblemTest.cpp -o BaseProblemTest.o

# Make SmallBiobjectiveSPProblem.o
SmallBiobjectiveSPProblem.o: SmallBiobjectiveSPProblem.cpp SmallBiobjectiveSPProblem.h ../PointAndSolution.h ../PointAndSolution.cpp ../BaseProblem.h ../BaseProblem.cpp ../Point.h
	$(CC) $(CPPFLAGS) -c SmallBiobjectiveSPProblem.cpp -o SmallBiobjectiveSPProblem.o

# Make SmallTripleobjectiveSPProblem.o
SmallTripleobjectiveSPProblem.o: SmallTripleobjectiveSPProblem.cpp SmallTripleobjectiveSPProblem.h ../PointAndSolution.h ../PointAndSolution.cpp ../BaseProblem.h ../BaseProblem.cpp ../Point.h
	$(CC) $(CPPFLAGS) -c SmallTripleobjectiveSPProblem.cpp -o SmallTripleobjectiveSPProblem.o

# Make NonOptimalStartingPointsProblem.o
NonOptimalStartingPointsProblem.o: NonOptimalStartingPointsProblem.cpp NonOptimalStartingPointsProblem.h ../PointAndSolution.h ../PointAndSolution.cpp ../BaseProblem.h ../BaseProblem.cpp ../Point.h
	$(CC) $(CPPFLAGS) -c NonOptimalStartingPointsProblem.cpp -o NonOptimalStartingPointsProblem.o


# Make Point.o
Point.o: ../Point.h ../Point.cpp ../DifferentDimensionsException.h ../NegativeApproximationRatioException.h ../NotPositivePointException.h
	$(CC) $(CPPFLAGS) -c ../Point.cpp -o Point.o

# Make Line2D.o
Line2D.o: ../Line2D.h ../Line2D.cpp ../Point.h Point.o ../Not2DPointsException.h ../SamePointsException.h ../VerticalLineException.h ../ParallelLinesException.h
	$(CC) $(CPPFLAGS) -c ../Line2D.cpp -o Line2D.o

# Make Hyperplane.o
Hyperplane.o: ../Hyperplane.h ../Hyperplane.cpp ../Point.h Point.o ../DifferentDimensionsException.h ../SamePointsException.h ../Not2DPointsException.h ../Not2DHyperplanesException.h ../NonExistentCoefficientException.h ../ParallelHyperplanesException.h
	$(CC) $(CPPFLAGS) -c ../Hyperplane.cpp -o Hyperplane.o

# Clean object files and executables
clean: 
	rm -f Point.o Line2D.o Hyperplane.o BaseProblemTest.o SmallBiobjectiveSPProblem.o SmallTripleobjectiveSPProblem.o NonOptimalStartingPointsProblem.o PointTest.out Line2DTest.out HyperplaneTest.out BaseProblemTest.out

