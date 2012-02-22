# Makefile for all unit tests.
# 
# Available unit tests are currently:
# - PointTest.cpp
# - Line2DTest.cpp
# - BaseProblemTest.cpp  &  VerySmallProblem.h  &  VerySmallProblem.cpp  &  
#   SmallGraphProblem.h  &  SmallGraphProblem.cpp
#


CC=g++
CPPFLAGS=-Wall
CPPLIBS=-lgtest


# Make all unit tests
all: PointTest.out Line2DTest.out BaseProblemTest.out

# Make PointTest.out
PointTest.out: PointTest.cpp ../Point.h Point.o ../DifferentDimensionsException.h
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o PointTest.cpp -o PointTest.out

# Make Line2DTest.out
Line2DTest.out: Line2DTest.cpp ../Point.h Point.o ../Line2D.h Line2D.o ../Not2DPointsException.h ../SamePointsException.h ../VerticalLineException.h ../ParallelLinesException.h
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o Line2D.o Line2DTest.cpp -o Line2DTest.out

# Make BaseProblemTest.out
BaseProblemTest.out: BaseProblemTest.cpp ../Point.h Point.o ../Line2D.h Line2D.o ../BaseProblem.h ../BaseProblem.cpp VerySmallProblem.h VerySmallProblem.cpp SmallGraphProblem.h SmallGraphProblem.cpp 
	$(CC) $(CPPFLAGS) $(CPPLIBS) Point.o Line2D.o VerySmallProblem.cpp SmallGraphProblem.cpp BaseProblemTest.cpp -o BaseProblemTest.out


# Make Point.o
Point.o: ../Point.h ../Point.cpp ../DifferentDimensionsException.h
	$(CC) $(CPPFLAGS) -c ../Point.cpp -o Point.o

# Make Line2D.o
Line2D.o: ../Line2D.h ../Line2D.cpp ../Point.h Point.o ../Not2DPointsException.h ../SamePointsException.h ../VerticalLineException.h ../ParallelLinesException.h
	$(CC) $(CPPFLAGS) -c ../Line2D.cpp -o Line2D.o

