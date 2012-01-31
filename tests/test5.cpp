/* test5.cpp */


#include <iostream>

#include "../chord.h"
#include "../ParallelLinesException.h"
#include "../SamePointsException.h"


using std::cout;
using std::endl;

using pareto_approximator::SamePointsException;
using pareto_approximator::ParallelLinesException;


int 
main(void)
{
  try {
    throw SamePointsException();
  }
  catch (SamePointsException e) {
    cout << e.what() << endl;
  }
  try {
    throw ParallelLinesException();
  }
  catch (ParallelLinesException e) {
    cout << e.what() << endl;
  }

  return 0;
}

