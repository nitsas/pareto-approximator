/*! \file utility.cpp
 *  \brief The implementation of some utility functions.
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  Won't `#include` "utility.h". In fact, "utility.h" will `#include` 
 *  "utility.cpp" because it contains function templates (which don't allow 
 *  us to split declaration from definition).
 */


#include <iostream>
#include <unistd.h>
#include <sys/wait.h>
#include <string>
#include <assert.h>

#include "NonDominatedSet.h"


/*!
 *  \weakgroup ParetoApproximator Everything needed for the chord algorithm.
 *  @{
 */


// An unnamed namespace containing helper functions. (declarations here)
// (implementations at the end of this file)
namespace {


using pareto_approximator::PointAndSolution;
using pareto_approximator::Facet;


//! Make qconvex's input file. The given points will be the input.
/*!
 *  \param points The points whose convex hull we need qconvex to compute
 *                for us.
 *  \param filename The name of the file we will create. (qconvex's input file)
 *  \param spaceDimension The dimension of the space that the points live in.
 *  
 *  \sa pareto_approximator::computeConvexHull()
 */
template <class S> 
void 
writePointsToQconvexInputFile(
                  std::vector< PointAndSolution<S> > points, 
                  std::string filename, unsigned int spaceDimension);


//! Parse qconvex's output file. Return a list of the convex hull's facets.
/*!
 *  \param filename The name of the file we will parse. (qconvex's output file)
 *  \param points The points whose convex hull we asked qconvex to compute
 *                for us. Their order counts - it must be the same as their 
 *                order in qconvex's input file.
 *  \param spaceDimension The dimension of the space that the points live in.
 *  
 *  \sa pareto_approximator::computeConvexHull()
 */
template <class S> 
std::list< Facet<S> > 
readFacetsFromQconvexOutputFile(std::string filename, 
                    std::vector< PointAndSolution<S> > points, 
                    unsigned int spaceDimension);


//! Normalizes a vector of double. (in place)
/*!
 *  \param v A vector (as a std::vector<double>).
 *
 *  After the operation the vector will be normalized, i.e. its magnitude 
 *  (in other words length or L2-norm) will be 1.
 */
void 
normalizeVector(std::vector<double> & v);


}  // namespace


//! The namespace containing everything needed for the chord algorithm.
namespace pareto_approximator {


//! Compute the convex hull of the given set of points.
/*!
 *  \param points A (const reference to a) std::vector of points. 
 *                (PointAndSolution<S> instances)
 *  
 *  \sa BaseProblem::doCraft()
 */
template <class S> 
std::list< Facet<S> > 
computeConvexHull(const std::vector< PointAndSolution<S> > & points, 
                  unsigned int spaceDimension)
{
  // The files we'll use to interface with qconvex.
  std::string qconvexInputFilename = "qconvex-input.txt";
  std::string qconvexOutputFilename = "qconvex-output.txt";

  // First make qconvex's input file:
  writePointsToQconvexInputFile(points, qconvexInputFilename, spaceDimension);

  std::list< Facet<S> > facets;

  // Make a subprocess (child) that will exec qconvex to compute 
  // the convex hull:
  pid_t pid = fork();
  if (pid < 0) {
    // could not fork()
    std::cout << "Failed to fork... Exiting" << std::endl;
    exit(-1);
  }
  else if (pid == 0) {
    // child process 
    // run qconvex (input: qconvex-input.txt, output: qconvex-output.txt)
    int rv = execlp("qconvex", "qconvex", "i", "n", "Qt", 
                    "PF0.000000000000000000000001", 
                    "TI", "qconvex-input.txt", 
                    "TO", "qconvex-output.txt", NULL);
    // Added option "PFn" where n is a lower bound for the area of 
    // facets to be printed (facets with area < n will not be printed) to 
    // avoid degenerate facets with zero area.
    // - e.g. "PF0.000000000000000000000001" will not print facets with 
    //   area less than 0.000000000000000000001
    // - degenerate facets with zero area might appear due to option "Qt"
    //   (how come? all the vertices might belong to the same ridge 
    //   of the original non simplicial facet)

    if (rv == -1) {
      std::cout << "An error occured while trying to call qconvex... Exiting" 
                << std::endl;
      exit(-1);
    }

    // Code will have ended by here. 
    // - Either via qconvex exiting without error or the exit(-1) in case a 
    //   qconvex error occurs. 
    // - Won't reach the return statement.
  }
  else {
    // parent process
    int childStatus;
    waitpid(pid, &childStatus, 0);
    if ( not WIFEXITED(childStatus) ) {
      if ( WIFSIGNALED(childStatus) ) {
        std::cout << "ERROR: qconvex got a signal that caused it to exit (" 
                  << WTERMSIG(childStatus) << "). Exiting" << std::endl;
        exit(-1);
      }
      else {
        std::cout << "ERROR: qconvex did not exit normally. Exiting"
                  << std::endl;
        exit(-1);
      }
    }
    else if ( WEXITSTATUS(childStatus) != 0 ) {
      std::cout << "ERROR: qconvex exited with errorcode: " 
                << WEXITSTATUS(childStatus) << std::endl
                << "Exiting" << std::endl;
      exit(-1);
    }
    else {
      // parse qconvex's output file (qconvex-output.txt) and make the facets
      facets = readFacetsFromQconvexOutputFile(qconvexOutputFilename, 
                                               points, 
                                               spaceDimension);
      // return the facets outside the ifs (to avoid compiler warnings 
      // about reaching the end of non-void function)
    }
  }

  // Only the parent will get here and only after succesfully reading 
  // the facets from qconvex's output file. (qconvex will have exited 
  // without error as well)
  return facets;
}


/*!
 *  \brief Filter a sequence of PointAndSolution instances and return 
 *         only the non-dominated ones.
 *
 *  \param first An iterator to the first element in the sequence.
 *  \param last An iterator to the past-the-end element in the sequence.
 *  
 *  We will use a pareto_approximator::NonDominatedSet to discard dominated 
 *  points.
 *  
 *  \sa NonDominatedSet
 */
template <class S> 
std::vector< PointAndSolution<S> > 
filterDominatedPoints(
      typename std::vector< PointAndSolution<S> >::const_iterator first, 
      typename std::vector< PointAndSolution<S> >::const_iterator last)
{
  NonDominatedSet< PointAndSolution<S> > filter(first, last);
  return std::vector< PointAndSolution<S> >(filter.begin(), filter.end());
}


//! Discard facets not useful for generating new Pareto points.
/*!
 *  \param facets A (reference to a) list of facets.
 *  
 *  Discard facets with all normal vector coefficients non-positive (<= 0).
 *  Facets with no positive normal vector coefficient are not useful for 
 *  generating new Pareto optimal points.
 *  
 *  Only facets with all-positive or mixed (i.e. containing at least some 
 *  positive coefficients) normal vectors can be used to generate new 
 *  Pareto optimal points.
 *  
 *  \sa Facet, BaseProblem::doChord() and BaseProblem::doCraft()
 */
template <class S> 
void
discardUselessFacets(std::list< Facet<S> > & facets)
{
  typename std::list< Facet<S> >::iterator it;

  for (it = facets.begin(); it != facets.end(); ) 
    if (it->hasAllNormalVectorElementsNonPositive())
      it = facets.erase(it);
    else
      ++it;
}


/*! \brief Choose the Facet instance with the largest local approximation 
 *         error upper bound from sequence of Facet instances.
 *  
 *  \param first An iterator to the first element in the sequence.
 *  \param last An iterator to the past-the-end element in the sequence.
 *  \return An iterator to the first element in the range that has the  
 *          largest local approximation error upper bound. If no element 
 *          is a non-boundary facet the function returns "last".
 *  
 *  The Facet::getLocalApproximationErrorUpperBound() method is used (of 
 *  course) for the facet's local approximation error upper bound.
 *  
 *  Boundary facets (i.e. those with isBoundaryFacet() == true) are ignored.
 *  
 *  If all the facets in the sequence are boundary facets the iterator 
 *  "last" is returned.
 *  
 *  \sa Facet
 */
template <class S> 
typename std::list< Facet<S> >::iterator 
chooseFacetWithLargestLocalApproximationErrorUpperBound(
                  typename std::list< Facet<S> >::iterator first, 
                  typename std::list< Facet<S> >::iterator last)
{
  typename std::list< Facet<S> >::iterator it, max = last;

  for (it = first; it != last; ++it) {
    // Is it a non-boundary facet?
    if (not it->isBoundaryFacet()) 
      // Is it the first non-boundary facet?
      // or 
      // Does it have a larger local approximation error upper bound 
      // than the one max has?
      if ( max == last or it->getLocalApproximationErrorUpperBound() > 
                          max->getLocalApproximationErrorUpperBound() ) 
        max = it;
      // else ignore it
    // else ignore it
  }

  return max;
}


//! Generate a weight vector (for comb()) using the given facet.
/*!
 *  \param facet A Facet instance. (Its vertices' weightsUsed attributes 
 *               will be needed if the facet's normal vector is not 
 *               all-positive.)
 *  \return A std::vector of objective weights (for BaseProblem::comb()).
 *  
 *  The resulting weights will be:
 *  - Either the facet's normal vector. 
 *    (if it has no negative elements)
 *  - Or the mean of the weights used to obtain the facet's vertices.
 *  
 *  \sa BaseProblem, BaseProblem::comb() and 
 *      BaseProblem::generateNewParetoPoint()
 */
template <class S> 
std::vector<double> 
generateNewWeightVector(const Facet<S> & facet)
{
  std::vector<double> weights;

  if (facet.hasAllNormalVectorElementsNonNegative()) {
    // Use the facet's normal vector (i.e. the facet's slope) as weights.
    weights = facet.getNormalVector();
  }
  else {
    // Use the mean of the facet's vertex weights (weightsUsed) as weights.
    // - "weights" will be a std::vector<double> of weights W_{i}, where:
    //   /f$ W_{i} = sum_{j=1}^{j=facet.spaceDimension()} ( w_{ij} ) /f$,
    //   where w_{ij} is the i'th of the weights used to obtain the j'th 
    //   vertex of "facet".
    weights = facet.computeMeanVertexWeights();
  }
  normalizeVector(weights);

  return weights;
}


}  // namespace pareto_approximator


// An unnamed namespace containing helper functions. (implementations here)
// (declarations at the beginning of this file)
namespace {


//! Make qconvex's input file. The given points will be the input.
/*!
 *  \param points The points whose convex hull we need qconvex to compute
 *                for us.
 *  \param filename The name of the file we will create. (qconvex's input file)
 *  \param spaceDimension The dimension of the space that the points live in.
 *  
 *  \sa pareto_approximator::computeConvexHull()
 */
template <class S> 
void 
writePointsToQconvexInputFile(
                  std::vector< PointAndSolution<S> > points, 
                  std::string filename, unsigned int spaceDimension)
{
  // qcif stands for "QConvex's Input File"
  std::ofstream qcif(filename.c_str(), std::ios::out | std::ios::trunc);

  if (not qcif.is_open()) {
    std::cout << "An error occured while opening file \"" << filename 
              << "\" for output... Exiting" << std::endl;
    exit(-1);
  }
  // else 

  // make qconvex's input file
  qcif << spaceDimension << std::endl;
  qcif << points.size() << std::endl;
  typename std::vector< PointAndSolution<S> >::iterator pit;
  for (pit = points.begin(); pit != points.end(); ++pit)
    qcif << pit->point << std::endl;

  qcif.close();
}


//! Parse qconvex's output file. Return a list of the convex hull's facets.
/*!
 *  \param filename The name of the file we will parse. (qconvex's output file)
 *  \param points The points whose convex hull we asked qconvex to compute
 *                for us. Their order counts - it must be the same as their 
 *                order in qconvex's input file.
 *  \param spaceDimension The dimension of the space that the points live in.
 *  
 *  \sa pareto_approximator::computeConvexHull()
 */
template <class S> 
std::list< Facet<S> > 
readFacetsFromQconvexOutputFile(std::string filename, 
                    std::vector< PointAndSolution<S> > points, 
                    unsigned int spaceDimension)
{
  std::list< Facet<S> > facets;

  // qcof stands for "QConvex's Output File"
  std::ifstream qcof(filename.c_str());

  if (not qcof.is_open()) {
    std::cout << "An error occured while opening file \"" << filename
              << "\" for input... Exiting" << std::endl;
    exit(-1);
  }
  // else

  // start parsing qconvex's output file
  unsigned int numFacets, i, j, vertexIndex;
  qcof >> numFacets;

  if (numFacets == 0)
    // nothing else to read
    return facets;
  
  // first read each facet's vertices
  std::vector< typename Facet<S>::VerticesVector > vertexSets(numFacets);
  for (i = 0; i != numFacets; ++i) {
    for (j = 0; j != spaceDimension; ++j) {
      qcof >> vertexIndex;
      // Make sure that the vertex index is a valid index.
      // - If the vertex index is greater than all valid indices it means 
      //   that qconvex encountered a non-simplicial facet and we had not 
      //   told qconvex to triangulate such facets.
      // - We call a facet in a d-dimensional space non-simplicial when it 
      //   consists of more than d vertices.
      // - Triangulating a non-simplicial facet with v vertices that lives 
      //   in a d-dimensional space (v > d) is the process of splitting the 
      //   (non-simplicial) facet into a set of coplanar simplicial facets.
      //   (each of the resulting simplicial facets will have d vertices and 
      //   they will all have the same normal vector)
      // - qconvex should triangulate non-simplicial facets if we give it 
      //   the "Qt" option (along with the "i" option) when calling it.
      // - Not sure (yet) if qconvex may encounter non-simplicial facets 
      //   in this use case i.e. convex hull of points belonging to a 
      //   Pareto set. (this only applies to the case when the given set of 
      //   points only contains Pareto optimal points)
      assert(vertexIndex < points.size());
      vertexSets[i].push_back(points[vertexIndex]);
    }
  }

  unsigned int temp;
  // consume (read) the next line 
  // - contains the facet normal vector size + 1
  qcof >> temp;
  assert(temp == spaceDimension + 1);

  unsigned int numOfNormalVectors;
  // consume (read) the next line
  // - contains the number of normal vectors (number of facets)
  qcof >> numOfNormalVectors;
  assert(numOfNormalVectors == numFacets);

  // now read each facet's normal vector and make the Facet instance
  // - qconvex facet's normal vectors face outwards (from the convex hull)
  //   but ours face inwards so we will reverse each normal vector's sign
  double normalVectorElement, facetOffset;
  for (i = 0; i != numFacets; ++i) {
    std::vector<double> facetNormal(spaceDimension);
    for (j = 0; j != spaceDimension; ++j) {
      qcof >> normalVectorElement;
      facetNormal[j] = - normalVectorElement;
    }
    // consume the next number (it's the facet offset - we don't need it)
    qcof >> facetOffset;
    // make the facet and push it onto the list of facets
    facets.push_back(Facet<S>(vertexSets[i].begin(), vertexSets[i].end(), 
                              facetNormal.begin(), facetNormal.end()));
  }

  // We should be at the end of the file. Make sure:
  double consume = 0.0;
  qcof >> consume;
  assert(qcof.eof());

  // Close the file and return the list of facets:
  qcof.close();

  return facets;
}


//! Normalizes a vector of double. (in place)
/*!
 *  \param v A vector (as a std::vector<double>).
 *
 *  After the operation the vector will be normalized, i.e. its magnitude 
 *  (in other words length or L2-norm) will be 1.
 */
void 
normalizeVector(std::vector<double> & v) 
{
  double l2Norm = arma::norm(arma::vec(v), 2);
  
  for (unsigned int i = 0; i != v.size(); ++i)
    v[i] /= l2Norm;
}


}  // namespace


/*! @} */
