/*! \file NonDominatedSetTest.cpp
 *  \brief Unit test for the NonDominatedSet class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <string>
#include <vector>

#include "gtest/gtest.h"
#include "../Point.h"
#include "../PointAndSolution.h"
#include "../NonDominatedSet.h"


using pareto_approximator::Point;
using pareto_approximator::PointAndSolution;
using pareto_approximator::NonDominatedSet;


namespace {


// The fixture for testing the NonDominatedSet class template.
class NonDominatedSetTest : public ::testing::Test 
{
  protected:
    NonDominatedSetTest() { }
    ~NonDominatedSetTest() { }

    virtual void SetUp()
    {
      std::vector< Point > elements;
      elements.push_back(Point(1, 2));
      elements.push_back(Point(2, 1));
      nds.insert(elements.begin(), elements.end());
    }

    NonDominatedSet<Point> emptyNds;
    NonDominatedSet<Point> nds;
};


// Test NonDominatedSet's constructors.
TEST_F(NonDominatedSetTest, NonDominatedSetBeginEndAndConstructorsWork)
{
  EXPECT_EQ(emptyNds.size(), 0);
  EXPECT_TRUE(emptyNds.begin() == emptyNds.end());
  
  EXPECT_EQ(nds.size(), 2);
  EXPECT_FALSE(nds.begin() == nds.end());
  NonDominatedSet<Point>::iterator it = nds.begin();
  EXPECT_EQ(*it, Point(1, 2));
  ++it;
  EXPECT_EQ(*it, Point(2, 1));
  ++it;
  EXPECT_TRUE(it == nds.end());
}


TEST_F(NonDominatedSetTest, NonDominatedSetSimpleMethodsWork)
{
  EXPECT_TRUE(emptyNds.empty());
  EXPECT_EQ(emptyNds.size(), 0);
  EXPECT_TRUE(emptyNds.find(Point(1, 1)) == emptyNds.end());
  EXPECT_FALSE(emptyNds.dominates(Point(1, 1)));
  emptyNds.clear();
  EXPECT_TRUE(emptyNds.empty());
  EXPECT_EQ(emptyNds.size(), 0);

  EXPECT_FALSE(nds.empty());
  EXPECT_EQ(nds.size(), 2);
  NonDominatedSet<Point>::iterator it = nds.begin();
  EXPECT_TRUE(nds.find(Point(1, 2)) == it);
  ++it;
  EXPECT_TRUE(nds.find(Point(2, 1)) == it);
  EXPECT_TRUE(nds.find(Point(3, 3)) == nds.end());
  EXPECT_TRUE(nds.dominates(Point(1, 4)));
  EXPECT_FALSE(nds.dominates(Point(1, 1)));
  nds.clear();
  EXPECT_TRUE(nds.empty());
  EXPECT_EQ(nds.size(), 0);
}


TEST_F(NonDominatedSetTest, NonDominatedSetInsertWorks)
{
  Point p1(1, 2);
  Point p2(2, 1);
  Point p3(3, 3);
  std::vector<Point> elements;
  elements.push_back(p1);
  elements.push_back(p2);
  elements.push_back(p3);
  emptyNds.insert(elements.begin(), elements.end());
  EXPECT_EQ(emptyNds.size(), 2);
  NonDominatedSet<Point>::iterator it = emptyNds.begin();
  EXPECT_EQ(*it, Point(1, 2));
  ++it;
  EXPECT_EQ(*it, Point(2, 1));
  ++it;
  EXPECT_TRUE(it == emptyNds.end());

  Point p4(1, 4);
  Point p5(1, 2);
  Point p6(1.5, 1.5);
  nds.insert(p4);
  EXPECT_EQ(nds.size(), 2);
  nds.insert(p5);
  EXPECT_EQ(nds.size(), 2);
  nds.insert(p6);
  EXPECT_EQ(nds.size(), 3);
  NonDominatedSet<Point>::iterator it2 = nds.begin();
  EXPECT_EQ(*it2, Point(1, 2));
  ++it2;
  EXPECT_EQ(*it2, Point(1.5, 1.5));
  ++it2;
  EXPECT_EQ(*it2, Point(2, 1));
  ++it2;
  EXPECT_TRUE(it2 == nds.end());
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

