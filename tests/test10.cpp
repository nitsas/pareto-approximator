/* test10.cpp */


#include <iostream>

#include "../Point.h"
#include "../Line2D.h"


using std::cout;
using std::endl;

using namespace pareto_approximator;


int 
main(void)
{
  Point p1(0, -2);
  Point p2(3, 1);
  Line2D l1(p1, p2);
  cout << l1 << endl;

  return 0;
}

