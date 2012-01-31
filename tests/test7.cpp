/* test7.cpp */


#include <iostream>
#include <sstream>
#include <string>

#include "../Point.h"


using std::cout;
using std::endl;
using std::stringstream;

int 
main(void)
{
  std::stringstream ss;

  ss << "Hello world! " << 557 << " lalala";
  const char* msg = ss.str().c_str();
  cout << msg << endl;

  double xx = 4.0;

  return 0;
}

