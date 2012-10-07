# Makefile for all unit tests.
# 
# Available unit tests are currently:
# - PointTest.cpp
# - HyperplaneTest.cpp
# - BaseProblemTest.cpp  &  NonOptimalStartingPointsProblem.h  &  
#		NonOptimalStartingPointsProblem.cpp  &  SmallBiobjectiveSPProblem.h  &  
#		SmallBiobjectiveSPProblem.cpp $ SmallTripleobjectiveSPProblem.h & 
#   SmallTrimpleobjectiveSPProblem.cpp & 
#   TripleobjectiveWithNegativeWeightsProblem.cpp & 
#   TripleobjectiveWithNegativeWeightsProblem.h
# 
# Author:  Christos Nitsas
# Date:    2012
#


CC=g++
CPPFLAGS=-Wall -Wextra -Werror -g
CPPLIBS=-lgtest -larmadillo


# Make all unit tests
all: PointTest.out HyperplaneTest.out NonDominatedSetTest.out BaseProblemTest.out

# Make PointTest.out
PointTest.out: PointTest.cpp ../Point.h Point.o ../DifferentDimensionsException.h ../NegativeApproximationRatioException.h ../NotPositivePointException.h ../NonExistentCoordinateException.h
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o PointTest.cpp -o $@

# Make HyperplaneTest.out
HyperplaneTest.out: HyperplaneTest.cpp ../Point.h Point.o ../Hyperplane.h Hyperplane.o ../NonExistentCoefficientException.h
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o Hyperplane.o HyperplaneTest.cpp -o $@

# Make NonDominatedSetTest.out
NonDominatedSetTest.out: NonDominatedSetTest.o Point.o
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o NonDominatedSetTest.o -o $@

# Make BaseProblemTest.out
BaseProblemTest.out: BaseProblemTest.o SmallBiobjectiveSPProblem.o SmallTripleobjectiveSPProblem.o NonOptimalStartingPointsProblem.o TripleobjectiveWithNegativeWeightsProblem.o Point.o Hyperplane.o
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o Hyperplane.o NonOptimalStartingPointsProblem.o SmallBiobjectiveSPProblem.o SmallTripleobjectiveSPProblem.o TripleobjectiveWithNegativeWeightsProblem.o BaseProblemTest.o -o $@


# Make BaseProblemTest.o
BaseProblemTest.o: BaseProblemTest.cpp ../Point.h ../PointAndSolution.h ../PointAndSolution.cpp NonOptimalStartingPointsProblem.h SmallBiobjectiveSPProblem.h SmallTripleobjectiveSPProblem.h TripleobjectiveWithNegativeWeightsProblem.h ../BaseProblem.h ../BaseProblem.cpp ../NonDominatedSet.h ../NonDominatedSet.cpp
	$(CC) $(CPPFLAGS) -c BaseProblemTest.cpp -o $@

# Make NonDominatedSetTest.o
NonDominatedSetTest.o: NonDominatedSetTest.cpp ../Point.h ../PointAndSolution.h ../PointAndSolution.cpp ../NonDominatedSet.h ../NonDominatedSet.cpp
	$(CC) $(CPPFLAGS) -c NonDominatedSetTest.cpp -o $@

# Make SmallBiobjectiveSPProblem.o
SmallBiobjectiveSPProblem.o: SmallBiobjectiveSPProblem.cpp SmallBiobjectiveSPProblem.h ../PointAndSolution.h ../PointAndSolution.cpp ../BaseProblem.h ../BaseProblem.cpp ../Point.h ../NonDominatedSet.h ../NonDominatedSet.cpp ../NotEnoughAnchorPointsException.h
	$(CC) $(CPPFLAGS) -c SmallBiobjectiveSPProblem.cpp -o $@

# Make SmallTripleobjectiveSPProblem.o
SmallTripleobjectiveSPProblem.o: SmallTripleobjectiveSPProblem.cpp SmallTripleobjectiveSPProblem.h ../PointAndSolution.h ../PointAndSolution.cpp ../BaseProblem.h ../BaseProblem.cpp ../Point.h ../NonDominatedSet.h ../NonDominatedSet.cpp ../NotEnoughAnchorPointsException.h
	$(CC) $(CPPFLAGS) -c SmallTripleobjectiveSPProblem.cpp -o $@

# Make NonOptimalStartingPointsProblem.o
NonOptimalStartingPointsProblem.o: NonOptimalStartingPointsProblem.cpp NonOptimalStartingPointsProblem.h ../PointAndSolution.h ../PointAndSolution.cpp ../BaseProblem.h ../BaseProblem.cpp ../Point.h ../NonDominatedSet.h ../NonDominatedSet.cpp ../NotEnoughAnchorPointsException.h
	$(CC) $(CPPFLAGS) -c NonOptimalStartingPointsProblem.cpp -o $@

# Make TripleobjectiveWithNegativeWeightsProblem.o
TripleobjectiveWithNegativeWeightsProblem.o: TripleobjectiveWithNegativeWeightsProblem.cpp TripleobjectiveWithNegativeWeightsProblem.h ../PointAndSolution.h ../PointAndSolution.cpp ../BaseProblem.h ../BaseProblem.cpp ../Point.h ../NonDominatedSet.h ../NonDominatedSet.cpp ../NotEnoughAnchorPointsException.h
	$(CC) $(CPPFLAGS) -c TripleobjectiveWithNegativeWeightsProblem.cpp -o $@

# Make Point.o
Point.o: ../Point.h ../Point.cpp ../DifferentDimensionsException.h ../NegativeApproximationRatioException.h ../NotPositivePointException.h
	$(CC) $(CPPFLAGS) -c ../Point.cpp -o $@

# Make Hyperplane.o
Hyperplane.o: ../Hyperplane.h ../Hyperplane.cpp ../Point.h Point.o ../DifferentDimensionsException.h ../SamePointsException.h ../NonExistentCoefficientException.h ../PointIsOnTheHyperplaneException.h
	$(CC) $(CPPFLAGS) -c ../Hyperplane.cpp -o $@

# Clean object files and executables
clean: 
	rm -f Point.o Hyperplane.o BaseProblemTest.o SmallBiobjectiveSPProblem.o SmallTripleobjectiveSPProblem.o NonOptimalStartingPointsProblem.o TripleobjectiveWithNegativeWeightsProblem.o NonDominatedSetTest.o PointTest.out HyperplaneTest.out BaseProblemTest.out NonDominatedSetTest.out

