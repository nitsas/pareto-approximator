/*! \file PointAndSolutionTest.cpp
 *  \brief Unit test for the PointAndSolution class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <string>
#include <vector>
#include <algorithm>

#include "gtest/gtest.h"
#include "../Point.h"
#include "../PointAndSolution.h"


using pareto_approximator::Point;
using pareto_approximator::PointAndSolution;
using pareto_approximator::NullObjectException;


namespace {


// The fixture for testing the PointAndSolution class template.
class PointAndSolutionTest : public ::testing::Test 
{
  protected:
    PointAndSolutionTest() 
    {
      nullPointAndSolution = PointAndSolution<std::string>();

      pas1 = PointAndSolution<std::string>(Point(2.0, 3.0, 4.0, 5.0), 
                                           "solution1");

      std::vector<double> weights(4, 0.5);
      pas2 = PointAndSolution<std::string>(Point(1.0, 2.0, 3.0, 4.0), 
                                           "solution2", 
                                           weights.begin(), 
                                           weights.end());
    }

    ~PointAndSolutionTest() { }

    PointAndSolution<std::string> nullPointAndSolution;
    PointAndSolution<std::string> pas1;
    PointAndSolution<std::string> pas2;
};


// Test PointAndSolution's constructors.
TEST_F(PointAndSolutionTest, PointAndSolutionConstructorsWork)
{
  // first for the nullPointAndSolution instance (null instance)
  EXPECT_TRUE(nullPointAndSolution.isNull());

  // next for the pas1 instance
  EXPECT_EQ(Point(2.0, 3.0, 4.0, 5.0), pas1.point);
  EXPECT_EQ(std::string("solution1"), pas1.solution);
  EXPECT_FALSE(pas1.isNull());
  EXPECT_TRUE(pas1.weightsUsed.empty());

  // and for the pas2 instance
  EXPECT_EQ(Point(1.0, 2.0, 3.0, 4.0), pas2.point);
  EXPECT_EQ(std::string("solution2"), pas2.solution);
  EXPECT_FALSE(pas2.isNull());
  std::vector<double> weights(4, 0.5);
  EXPECT_TRUE(std::equal(pas2.weightsUsed.begin(), 
                         pas2.weightsUsed.end(), 
                         weights.begin()));
}


// Test PointAndSolution's methods. (most are shortcuts to Point's methods)
TEST_F(PointAndSolutionTest, PointAndSolutionMethodsWork)
{
  // first for the nullPointAndSolution instance (null instance)
  EXPECT_TRUE(nullPointAndSolution.isNull());
  EXPECT_THROW(nullPointAndSolution.dominates(pas1), NullObjectException);
  EXPECT_THROW(nullPointAndSolution.dimension(), NullObjectException);

  // next for the pas1 instance
  EXPECT_FALSE(pas1.isNull());
  EXPECT_EQ(pas1.point.dominates(pas2.point), pas1.dominates(pas2));
  EXPECT_THROW(pas1.dominates(nullPointAndSolution), NullObjectException);
  EXPECT_EQ(pas1.point.dimension(), pas1.dimension());

  // and for the pas2 instance
  EXPECT_FALSE(pas2.isNull());
  EXPECT_EQ(pas2.point.dominates(pas1.point), pas2.dominates(pas1));
  EXPECT_THROW(pas2.dominates(nullPointAndSolution), NullObjectException);
  EXPECT_EQ(pas2.point.dimension(), pas2.dimension());
}


TEST_F(PointAndSolutionTest, PointAndSolutionOperatorsWork)
{
  // first for the nullPointAndSolution instance (null instance)
  EXPECT_TRUE(nullPointAndSolution == PointAndSolution<std::string>());
  EXPECT_FALSE(nullPointAndSolution == pas1);
  EXPECT_FALSE(nullPointAndSolution == pas2);
  EXPECT_THROW(nullPointAndSolution < pas1, NullObjectException);

  // next for the pas1 and pas2 instances
  PointAndSolution<std::string> pas3(Point(2.0, 3.0, 4.0, 5.0), "solution3");
  EXPECT_EQ(pas1.point == pas3.point, pas1 == pas3);
  EXPECT_EQ(pas1.point == pas2.point, pas1 == pas2);
  EXPECT_EQ(pas1.point < pas2.point, pas1 < pas2);
  EXPECT_EQ(pas2.point < pas1.point, pas2 < pas1);
  EXPECT_THROW(pas1 < nullPointAndSolution, NullObjectException);
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv) 
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

