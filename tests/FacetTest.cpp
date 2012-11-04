/*! \file FacetTest.cpp
 *  \brief Unit test for the Facet class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <vector>
#include <algorithm>
#include <string>

#include "gtest/gtest.h"
#include "../Facet.h"
#include "../Point.h"
#include "../PointAndSolution.h"
#include "../BoundaryFacetException.h"


using namespace pareto_approximator;


namespace {


// The fixture for testing class Facet.
class FacetTest : public ::testing::Test 
{
  protected:
    FacetTest()
    {
      // For the purpose of this test we will create two Facet instances, 
      // one that is not a boundary facet (has a unique Lower Distal Point 
      // (LDP)) and one that is (no unique Lower Distal Point). 

      regularFacet = makeRegularFacet();
      boundaryFacet = makeBoundaryFacet();
      verticesOnlyFacet = makeFacetUsingVerticesOnly();
    }

    ~FacetTest() { }

    Facet<std::string> makeRegularFacet() 
    {
      // regular facet == has a unique Lower Distal Point

      // We will create the facet's vertices and normal-vector manually 
      // before calling the constructor.

      // first make the facet's vertices
      PointAndSolution<std::string> v122(Point(1, 2, 2), "solution1");
      v122.weightsUsed.push_back(1.0);
      v122.weightsUsed.push_back(0.0);
      v122.weightsUsed.push_back(0.0);
      PointAndSolution<std::string> v212(Point(2, 1, 2), "solution2");
      v212.weightsUsed.push_back(0.0);
      v212.weightsUsed.push_back(1.0);
      v212.weightsUsed.push_back(0.0);
      PointAndSolution<std::string> v221(Point(2, 2, 1), "solution3");
      v221.weightsUsed.push_back(0.0);
      v221.weightsUsed.push_back(0.0);
      v221.weightsUsed.push_back(1.0);
      // We initialized the vertices' coordinates and weightsUsed attributes 
      // in such a way that the facet's Lower Distal Point will be (1, 1, 1).
      // i.e. The hyperplanes through the facet's vertices with normal vectors
      //      the respective weightsUsed vectors intersect at (1, 1, 1).

      std::vector< PointAndSolution<std::string> > vertices;
      vertices.push_back(v122);
      vertices.push_back(v212);
      vertices.push_back(v221);

      // then make the facet's normal vector
      std::vector<double> facetNormal(3, 1.0);

      // now make and return the facet
      Facet<std::string> regularFacet(vertices.begin(), vertices.end(), 
                                      facetNormal.begin(), facetNormal.end());

      return regularFacet;
    }

    Facet<std::string> makeBoundaryFacet() 
    {
      // boundary facet == no unique Lower Distal Point

      // We will create the facet's vertices and normal-vector manually 
      // before calling the constructor.

      // first make the facet's vertices
      PointAndSolution<std::string> v122(Point(1, 2, 2), "solution1");
      v122.weightsUsed.push_back(1.0);
      v122.weightsUsed.push_back(0.0);
      v122.weightsUsed.push_back(0.0);
      PointAndSolution<std::string> v313(Point(3, 1, 3), "solution2");
      v313.weightsUsed.push_back(0.0);
      v313.weightsUsed.push_back(1.0);
      v313.weightsUsed.push_back(0.0);
      PointAndSolution<std::string> v214(Point(2, 1, 4), "solution3");
      v214.weightsUsed.push_back(0.0);
      v214.weightsUsed.push_back(1.0);
      v214.weightsUsed.push_back(0.0);
      // We initialized the vertices' coordinates and weightsUsed attributes 
      // in such a way that the facet has no unique LowerDistalPoint.
      // i.e. The hyperplanes through the facet's vertices with normal vectors
      //      the respective weightsUsed vectors intersect at an infinite 
      //      number of points (the line {x=1, y=1, z=anything})

      std::vector< PointAndSolution<std::string> > vertices;
      vertices.push_back(v122);
      vertices.push_back(v313);
      vertices.push_back(v214);

      // then make the facet's normal vector
      std::vector<double> facetNormal(3);
      facetNormal[0] = -1.0;
      facetNormal[1] = -3.0;
      facetNormal[2] = -1.0;

      // now make and return the facet
      Facet<std::string> boundaryFacet(vertices.begin(), vertices.end(), 
                                       facetNormal.begin(), facetNormal.end());

      return boundaryFacet;
    }

    Facet<std::string> makeFacetUsingVerticesOnly() 
    {
      // It will be a boundary facet as well.

      // boundary facet == no unique Lower Distal Point

      // We will only create the facet's vertices manually and then the 
      // constructor will compute and set the facet's normal vector.

      // first make the facet's vertices
      PointAndSolution<std::string> v12(Point(1, 2), "solution1");
      v12.weightsUsed.push_back(1.0);
      v12.weightsUsed.push_back(0.0);
      PointAndSolution<std::string> v11(Point(1, 1), "solution2");
      v11.weightsUsed.push_back(1.0);
      v11.weightsUsed.push_back(0.0);
      // We initialized the vertices' coordinates and weightsUsed attributes 
      // in such a way that the facet has no unique LowerDistalPoint.
      // i.e. The hyperplanes (in this case lines) through the facet's 
      //      vertices with normal vectors the respective weightsUsed vectors 
      //      intersect at an infinite number of points (the vertical line 
      //      {x=1, y=anything})

      std::vector< PointAndSolution<std::string> > vertices;
      vertices.push_back(v12);
      vertices.push_back(v11);

      // now call the constructor to make the facet
      Facet<std::string> verticesOnlyFacet(vertices.begin(), vertices.end());

      return verticesOnlyFacet;
    }

    // A regular facet. Has a unique Lower Distal Point (LDP).
    Facet<std::string> regularFacet;
    // A boundary facet. Hasn't got a unique LDP.
    Facet<std::string> boundaryFacet;
    // A (boundary) facet constructed using the vertices-only constructor.
    // The constructor computed the facet's normal. Hasn't got a unique LDP.
    Facet<std::string> verticesOnlyFacet;
};


// Test that Facet's constructor and accessors work as expected.
TEST_F(FacetTest, FacetConstructorsAndAccessorsWork) 
{
  // The constructor was called inside FacetTest's SetUp method.

  // first for the regularFacet instance
  EXPECT_EQ(3, regularFacet.spaceDimension());
  EXPECT_FALSE(regularFacet.isBoundaryFacet());
  EXPECT_EQ(2.0 / 3.0, regularFacet.getLocalApproximationErrorUpperBound());

  // now for the boundaryFacet instance
  EXPECT_EQ(3, boundaryFacet.spaceDimension());
  EXPECT_TRUE(boundaryFacet.isBoundaryFacet());
  EXPECT_THROW(boundaryFacet.getLocalApproximationErrorUpperBound(), 
               BoundaryFacetException);

  // now check verticesOnlyFacet's attributes
  EXPECT_EQ(2, verticesOnlyFacet.spaceDimension());
  EXPECT_TRUE(verticesOnlyFacet.isBoundaryFacet());
  EXPECT_THROW(verticesOnlyFacet.getLocalApproximationErrorUpperBound(), 
               BoundaryFacetException);
}


// Test that Facet's iterators and methods returning iterators work.
TEST_F(FacetTest, FacetIteratorsWork)
{
  // first for the regularFacet
  EXPECT_EQ(3.0, std::distance(regularFacet.beginVertex(), 
                               regularFacet.endVertex()));

  Facet<std::string>::ConstVertexIterator cvi = regularFacet.beginVertex();
  EXPECT_EQ(Point(1, 2, 2), cvi->point);
  EXPECT_EQ("solution1", cvi->solution);
  EXPECT_EQ(3, cvi->weightsUsed.size());
  EXPECT_EQ(1.0, cvi->weightsUsed[0]);
  EXPECT_EQ(0.0, cvi->weightsUsed[1]);
  EXPECT_EQ(0.0, cvi->weightsUsed[2]);
  ++cvi;
  EXPECT_EQ(Point(2, 1, 2), cvi->point);
  EXPECT_EQ("solution2", cvi->solution);
  EXPECT_EQ(3, cvi->weightsUsed.size());
  EXPECT_EQ(0.0, cvi->weightsUsed[0]);
  EXPECT_EQ(1.0, cvi->weightsUsed[1]);
  EXPECT_EQ(0.0, cvi->weightsUsed[2]);
  ++cvi;
  EXPECT_EQ(Point(2, 2, 1), cvi->point);
  EXPECT_EQ("solution3", cvi->solution);
  EXPECT_EQ(3, cvi->weightsUsed.size());
  EXPECT_EQ(0.0, cvi->weightsUsed[0]);
  EXPECT_EQ(0.0, cvi->weightsUsed[1]);
  EXPECT_EQ(1.0, cvi->weightsUsed[2]);
  ++cvi;
  EXPECT_TRUE(cvi == regularFacet.endVertex());

  // now for the verticesOnlyFacet
  EXPECT_EQ(2.0, std::distance(verticesOnlyFacet.beginVertex(), 
                               verticesOnlyFacet.endVertex()));
  cvi = verticesOnlyFacet.beginVertex();
  EXPECT_EQ(Point(1, 2), cvi->point);
  EXPECT_EQ("solution1", cvi->solution);
  EXPECT_EQ(2, cvi->weightsUsed.size());
  EXPECT_EQ(1.0, cvi->weightsUsed[0]);
  EXPECT_EQ(0.0, cvi->weightsUsed[1]);
  ++cvi;
  EXPECT_EQ(Point(1, 1), cvi->point);
  EXPECT_EQ("solution2", cvi->solution);
  EXPECT_EQ(2, cvi->weightsUsed.size());
  EXPECT_EQ(1.0, cvi->weightsUsed[0]);
  EXPECT_EQ(0.0, cvi->weightsUsed[1]);
  ++cvi;
  EXPECT_TRUE(cvi == verticesOnlyFacet.endVertex());
}


// Test that Facet::computeMeanVertexWeights() works.
TEST_F(FacetTest, ComputeMeanVertexWeightsWorks)
{
  // first for the regularFacet
  std::vector<double> rfmvw = regularFacet.computeMeanVertexWeights();
  EXPECT_EQ(3, rfmvw.size());
  EXPECT_EQ(1.0 / 3.0, rfmvw[0]);
  EXPECT_EQ(1.0 / 3.0, rfmvw[1]);
  EXPECT_EQ(1.0 / 3.0, rfmvw[2]);

  // now for the boundaryFacet
  std::vector<double> bfmvw = boundaryFacet.computeMeanVertexWeights();
  EXPECT_EQ(3, bfmvw.size());
  EXPECT_EQ(1.0 / 3.0, bfmvw[0]);
  EXPECT_EQ(2.0 / 3.0, bfmvw[1]);
  EXPECT_EQ(0.0, bfmvw[2]);
}


// Test that Facet::computeLowerDistalPoint() works.
TEST_F(FacetTest, ComputeLowerDistalPointWorks)
{
  // We implicitly tested it inside the FacetConstructorsAndAccessorsWork 
  // test as well. (Facet's constructor uses the aforementioned method to 
  // initialize some attributes)

  ASSERT_TRUE(regularFacet.computeLowerDistalPoint() != NULL);
  EXPECT_EQ(Point(1, 1, 1), *regularFacet.computeLowerDistalPoint());

  EXPECT_TRUE(boundaryFacet.computeLowerDistalPoint() == NULL);
}


// Test that Facet::ratioDistance() works.
TEST_F(FacetTest, RatioDistanceWorks)
{
  // We implicitly tested it inside the FacetConstructorsAndAccessorsWork 
  // test as well. (Facet's constructor uses the aforementioned method to 
  // initialize some attributes)

  EXPECT_EQ(2.0 / 3.0, regularFacet.ratioDistance(Point(1, 1, 1)));
  EXPECT_EQ(1.0 / 4.0, regularFacet.ratioDistance(Point(1, 2, 1)));
  EXPECT_EQ(0.0, regularFacet.ratioDistance(Point(8, 8, 8)));

  EXPECT_EQ(4.0 / 5.0, boundaryFacet.ratioDistance(Point(1, 1, 1)));
  EXPECT_EQ(3.0 / 6.0, boundaryFacet.ratioDistance(Point(2, 1, 1)));
  EXPECT_EQ(0.0, boundaryFacet.ratioDistance(Point(8, 8, 8)));
}


// Test that Facet::hasAllNormalVectorElementsNonPositive() works.
TEST_F(FacetTest, HasAllNormalVectorElementsNonPositiveWorks)
{
  EXPECT_FALSE(regularFacet.hasAllNormalVectorElementsNonPositive());
  EXPECT_TRUE(boundaryFacet.hasAllNormalVectorElementsNonPositive());

  // verticesOnlyFacet's normal vector could be either all-elements-positive 
  // or all-elements-non-positive but we expect the constructor to choose 
  // the all-elements-positive one 
  // 
  // reminder: for each set of n points there are two possible n-hyperplanes
  //           that pass through all the points; the two hyperplanes' (and 
  //           the corresponding facets') normal vectors are opposite
  EXPECT_FALSE(verticesOnlyFacet.hasAllNormalVectorElementsNonPositive());
}


}  // namespace


// Run all tests
int 
main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

