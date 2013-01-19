/*! \file examples/tripleobjective_shortest_path/RandomGraphProblem.h
 *  \brief The declaration of the RandomGraphProblem class. 
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef EXAMPLE_CLASS_RANDOM_GRAPH_PROBLEM_H
#define EXAMPLE_CLASS_RANDOM_GRAPH_PROBLEM_H


#include <boost/graph/adjacency_list.hpp>

#include "tripleobjective_shortest_path_example_common.h"
#include "../../Point.h"
#include "../../PointAndSolution.h"
#include "../../BaseProblem.h"
#include "../../NonDominatedSet.h"


using pareto_approximator::Point;
using pareto_approximator::PointAndSolution;
using pareto_approximator::BaseProblem;
using pareto_approximator::NonDominatedSet;


/*!
 *  \addtogroup TripleobjectiveShortestPathExample An example tripleobjective shortest path problem.
 *  
 *  @{
 */


//! Everything needed for the example tripleobjective shortest path problem.
namespace tripleobjective_shortest_path_example {


//! A class representing a tripleobjective shortest path problem.
/*!
 *  Problem info
 *  --------------------
 *  A class representing a tripleobjective shortest path problem on an 
 *  undirected, random and with no parallel edges boost graph.
 *
 *  More info:
 *  - Edges have three (different) weights called "black", "red" and "green". 
 *  - Two vertices, the source (s) and target (t) are singled out. 
 *  - A problem solution is a path from s to t. 
 *
 *  The goal is to find a path that minimizes three objective functions, 
 *  which we will call Black, Red and Green. They are:
 *  - Black: The sum of "black" weights for all the path's edges.
 *  - Red: The sum of "red" weights for all the path's edges.
 *  - Green: The sum of "green" weights for all the path's edges.
 *  
 *  It is easy to see that every problem solution corresponds to a point 
 *  in objective space. RandomGraphProblem::computeConvexParetoSet() will 
 *  try to find a set of points whose convex combinations approximately 
 *  dominate every point in the problem's Pareto set (the set of 
 *  non-dominated points).
 *  
 *  Reminder: The user can set the degree of approximation when he calls 
 *  RandomGraphProblem::computeConvexParetoSet().
 *  
 *  What we had to do
 *  --------------------
 *  RandomGraphProblem instances inherit BaseProblem::computeConvexParetoSet() 
 *  directly from BaseProblem so the only thing we had to implement was 
 *  the RandomGraphProblem::comb() method (!) which is declared virtual in 
 *  BaseProblem. 
 *  
 *  We didn't have to implement anything else (except constructor/destructor). 
 *  RandomGraphProblem::makeGraph(), RandomGraphProblem::printPath() and 
 *  RandomGraphProblem::isTargetReachable() are just helpful, 
 *  problem-specific methods.
 *  
 *  /sa BaseProblem, PointAndSolution and Point
 */
class RandomGraphProblem : public BaseProblem<PredecessorMap>
{
  public:
    //! Constructor. Make a tripleobjective shortest path problem instance.
    /*!
     *  \param numVertices The number of vertices.
     *  \param numEdges The number of edges.
     *  \param minBlackWeight The lowest possible "black" edge weight.
     *  \param maxBlackWeight The maximum possible "black" edge weight.
     *  \param minRedWeight The lowest possible "red" edge weight.
     *  \param maxRedWeight The maximum possible "red" edge weight.
     *  \param minGreenWeight The lowest possible "green" edge weight.
     *  \param maxGreenWeight The maximum possible "green" edge weight.
     *  \param seed The random number generator's seed.
     *  
     *  A simple constructor for tripleobjective shortest path problems. (of the 
     *  type we described in RandomGraphProblem)
     *  - Makes an undirected random boost graph with no parallel edges. Three 
     *    integer weights on each edge called "black", "red" and "green".
     *  - "black" weights are random integers chosen uniformly from 
     *    [minBlackWeight, maxBlackWeight].
     *  - "red" weights are random integers chosen uniformly from 
     *    [minRedWeight, maxRedWeight].
     *  - "green" weights are random integers chosen uniformly from 
     *    [minGreenWeight, maxGreenWeight].
     *  - Three objective functions to minimize:
     *    + Black: the sum of all "black" weights in an s-t path.
     *    + Red: the sum of all "red" weights in an s-t path.
     *    + Green: the sum of all "green" weights in an s-t path.
     *  
     *  \sa ~RandomGraphProblem() and makeGraph()
     */
    RandomGraphProblem(int numVertices=100, int numEdges=800, 
                       int minBlackWeight=1, int maxBlackWeight=100, 
                       int minRedWeight=1, int maxRedWeight=100,
                       int minGreenWeight=1, int maxGreenWeight=100, 
                       int seed=1);

    //! Empty destructor.
    /*!
     *  \sa RandomGraphProblem()
     */
    ~RandomGraphProblem();

    //! Make an undirected random boost graph with no parallel edges. 
    /*!
     *  \param numVertices The number of vertices.
     *  \param numEdges The number of edges.
     *  \param minBlackWeight The lowest possible "black" edge weight.
     *  \param maxBlackWeight The maximum possible "black" edge weight.
     *  \param minRedWeight The lowest possible "red" edge weight.
     *  \param maxRedWeight The maximum possible "red" edge weight.
     *  \param minGreenWeight The lowest possible "green" edge weight.
     *  \param maxGreenWeight The maximum possible "green" edge weight.
     *  \param seed The random number generator's seed.
     *  
     *  - Three integer weights on each edge, called "black", "red" 
     *    and "green".
     *  - "black" weights are random integers chosen uniformly from 
     *    [minBlackWeight, maxBlackWeight].
     *  - "red" weights are random integers chosen uniformly from 
     *    [minRedWeight, maxRedWeight].
     *  - "green" weights are random integers chosen uniformly from 
     *    [minGreenWeight, maxGreenWeight].
     *  
     *  \sa RandomGraphProblem and RandomGraphProblem()
     */
    void makeGraph(int numVertices, int numEdges, 
                   int minBlackWeight, int maxBlackWeight, 
                   int minRedWeight, int maxRedWeight, 
                   int minGreenWeight, int maxGreenWeight, 
                   int seed);

    //! The comb routine we had to implement. 
    /*!
     *  \param first Iterator to the initial position in an 
     *               std::vector<double> containing the weights w_{i} of the 
     *               objectives (in the linear combination of objective 
     *               functions).
     *  \param last Iterator to the past-the-end position in an 
     *              std::vector<double> containing the weights w_{i} of the 
     *              objectives (in the linear combination of objective 
     *              functions).
     *  \return A PointAndSolution object containing an s-t path (P) that 
     *          minimizes \$f w_{0} * Black(P) + w_{1} * Red(P) + 
     *          w_{2} * Green(P) \$f and the corresponding point in objective 
     *          space.
     *  
     *  The vector of weights will only contain three weights for this 
     *  example, w_{0} (the Black objective function weight), w_{1} (the 
     *  Red objective function weight) and w_{2} (the Green objective 
     *  function weight.
     *  
     *  Minimizes linear combinations of the objective functions of the 
     *  following form:
     *  \$f w_{0} * Black(P) + w_{1} * Red(P) + w_{2} * Green(P) \$f,
     *  where P is an s-t path.
     *  
     *  \sa RandomGraphProblem and RandomGraphProblem::RandomGraphProblem().
     */
    PointAndSolution<PredecessorMap> comb(
                          std::vector<double>::const_iterator first, 
                          std::vector<double>::const_iterator last) const;

    //! Check if the target (t) is reachable.
    /*!
     *  \return True iff there is at least one path that connects source (s) 
     *          and target (t).
     *  \return False iff there is no path that connects source (s) and 
     *          target (t).
     */
    bool isTargetReachable();

    //! A simple method that prints s-t paths.
    /*!
     *  \param pred A map from each vertex to its predecessor in the path. 
     *              (Using a std::vector for the example.)
     *  
     *  Print the path as a series of vertex descriptors and edge weights 
     *  (formatted nicely).
     */
    void printPath(const PredecessorMap& pred) const;

    //! Find the exact Pareto set.
    /*!
     *  Will use boost's breadth_first_search and a custom visitor to 
     *  make something like the flood algorithm.
     */
    NonDominatedSet<Point> findExactParetoSet();

    //! Return a reference to the underlying graph.
    Graph& graph();
    //! Return a reference to the source vertex (s).
    Vertex& source();
    //! Return a reference to the target vertex (t).
    Vertex& target();
    //! Print the graph to a dot (Graphviz) file.
    void printGraphToDotFile(const char* filename="graph.dot");

  private:
    //! The underlying graph.
    Graph g_;
    //! The source vertex (s).
    Vertex s_;
    //! The target vertex (t).
    Vertex t_;
};


}  // namespace tripleobjective_shortest_path_example


/*! 
 *  @}
 */


#endif  // EXAMPLE_CLASS_RANDOM_GRAPH_PROBLEM_H
