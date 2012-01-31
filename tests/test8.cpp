/* test8.cpp */


#include <iostream>

#include "../DifferentDimensionsException.h"
#include "../Point.h"


using std::cout;
using std::cin;
using std::endl;

using namespace pareto_approximator;


int 
main(void)
{
  Point p1(4.5, 3.0, -2.0);
  Point p2(2.1, 5.72);

  try {
    throw DifferentDimensionsException(p1.dimension(), p2.dimension());
  }
  catch (DifferentDimensionsException& e) {
    cout << e.what() << endl;
  }

  Point p3;
  cin >> p3;
  cout << endl << "dimension: " << p3.dimension() << endl;
  cout << p3 << endl;

  return 0;
}

