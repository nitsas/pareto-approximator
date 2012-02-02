/* test9.cpp */


#include <iostream>

#include "../Point.h"


using std::cout;
using std::endl;
using std::ostream;

using namespace pareto_approximator;


class SimpleLine 
{
  public:
    SimpleLine() : m_(1), c_(0) { }
    SimpleLine(double m, double c) : m_(m), c_(c) { }

    static SimpleLine 
    lineThrough(Point& p1, Point& p2)
    {
      SimpleLine line( (p2.y - p1.y) / (p2.x - p1.x), p1.y - (p2.y - p1.y) * p1.x / (p2.x - p1.x) );
      return line;
    }

    friend ostream& operator<< (ostream& ostr, const SimpleLine& line)
    {
      return ostr << "y = " << line.m_ << " * x + (" << line.c_ << ")";
    }

    double m_;
    double c_;
}; 


int 
main(void)
{
  SimpleLine l1(0, 5);
  SimpleLine l2(1, 0);

  cout << l1 << endl << l2 << endl;

  Point p1(0, 0);
  Point p2(1, 1);

  SimpleLine l3 = SimpleLine::lineThrough(p1, p2);

  cout << endl << "Points:   " << p1 << "   and   " << p2 << endl;;
  cout << "Line through them:   " << l3 << endl;

  return 0;
}

