/*! \file RandomGraphProblem.h
 *  \brief The declaration of the RandomGraphProblem class. 
 */


#ifndef EXAMPLE_CLASS_RANDOM_GRAPH_PROBLEM_H
#define EXAMPLE_CLASS_RANDOM_GRAPH_PROBLEM_H


#include <string>
#include <boost/graph/adjacency_list.hpp>

#include "../../PointAndSolution.h"
#include "../../BaseProblem.h"


using std::string;

using pareto_approximator::PointAndSolution;
using pareto_approximator::BaseProblem;


//! A boost graph edge property. (bundled property)
class EdgeProperty 
{
  public:
    //! The edge's "black" weight.
    double black;
    //! The edge's "red" weight.
    double red;
};


//! A boost adjacency_list representing an undirected graph.
/*!
 *  - vertices are stored in a std::vector
 *  - no parallel edges
 *  - vertices have a color property (for BFS search)
 *  - edges have a bundled property (class EdgeProperty)
 */
typedef boost::adjacency_list<boost::setS, 
                              boost::vecS, 
                              boost::undirectedS, 
                              boost::no_property,
                              EdgeProperty>           Graph;
//! A graph vertex. 
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
//! A graph vertex iterator.
typedef boost::graph_traits<Graph>::vertex_iterator   VertexIterator;
//! A graph edge.
typedef boost::graph_traits<Graph>::edge_descriptor   Edge;
//! A graph edge iterator.
typedef boost::graph_traits<Graph>::edge_iterator     EdgeIterator;
//! A simple predecessor map.
/*!
 *  Maps a vertex to its predecessor. (e.g. in a shortest path)
 *  - Implemented using a std::vector of vertices since vertex descriptors 
 *    are integers. If the vertex storage type changes (in the Graph typedef) 
 *    this must change too.
 */
typedef std::vector<Vertex>                    PredecessorMap;


//! A class representing a biobjective shortest path problem.
/*!
 *  Problem info
 *  --------------------
 *  A class representing a biobjective shortest path problem on an 
 *  undirected, random and with no parallel edges boost graph.
 *
 *  More info:
 *  - Edges have two (different) weights called "black" and "red". 
 *  - Two vertices, the source (s) and target (t) are singled out. 
 *  - A problem solution is a path from s to t. 
 *
 *  The goal is to find a path that minimizes two objective functions, 
 *  which we will call Black and Red. They are:
 *  - Black: The sum of "black" weights for all the path's edges.
 *  - Red: The sum of "red" weights for all the path's edges.
 *  
 *  It is easy to see that every problem solution corresponds to a point 
 *  in objective space. RandomGraphProblem::approximateParetoSet() will 
 *  try to find a set of points whose convex combinations approximately 
 *  dominate every point in the problem's Pareto set (the set of 
 *  non-dominated points).
 *  
 *  Reminder: The user can set the degree of approximation when he calls 
 *  RandomGraphProblem::approximateParetoSet().
 *  
 *  What we had to do
 *  --------------------
 *  RandomGraphProblem instances inherit BaseProblem::approximateParetoSet() 
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
    //! Constructor. Make a biobjective shortest path problem instance.
    /*!
     *  \param numVertices The number of vertices.
     *  \param numEdges The number of edges.
     *  \param minBlackWeight The lowest possible "black" edge weight.
     *  \param maxBlackWeight The maximum possible "black" edge weight.
     *  \param minRedWeight The lowest possible "red" edge weight.
     *  \param maxRedWeight The maximum possible "red" edge weight.
     *  
     *  A simple constructor for biobjective shortest path problems. (of the 
     *  type we described in RandomGraphProblem)
     *  - Makes an undirected random boost graph with no parallel edges. Two 
     *    integer weights on each edge called "black" and "red".
     *  - "black" weights are random integers chosen uniformly from 
     *    [minBlackWeight, maxBlackWeight].
     *  - "red" weights are random integers chosen uniformly from 
     *    [minRedWeight, maxRedWeight].
     *  - Singles out two vertices (the first and last one created), the 
     *    source (s) and target (t).
     *  - Two objective functions to minimize:
     *    + Black: the sum of all "black" weights in an s-t path.
     *    + Red: the sum of all "red" weights in an s-t path.
     *  
     *  \sa ~RandomGraphProblem() and makeGraph()
     */
    RandomGraphProblem(int numVertices=1000, int numEdges=100000, 
                       int minBlackWeight=1, int maxBlackWeight=100, 
                       int minRedWeight=1, int maxRedWeight=100);

    //! Empty destructor.
    /*!
     *  \sa RandomGraphProblem()
     */
    ~RandomGraphProblem();

    //! The comb routine we had to implement. 
    /*!
     *  \param xWeight The x objective's (Black objective function) weight.
     *                 (in the linear combination of objective functions)
     *  \param yWeight The y objective's (Red objective function) weight.
     *                 (in the linear combination of objective functions)
     *  \return An s-t path (P) that minimizes \$f xWeight * Black(P) + 
     *          yWeight * Red(P) \$f and the corresponding point in 
     *          objective space.
     *  
     *  Minimizes linear combinations of the objective functions.
     *  
     *  \sa RandomGraphProblem, RandomGraphProblem() and BaseProblem::comb()
     */
    PointAndSolution<PredecessorMap> comb(double xWeight, double yWeight);

    //! Make an undirected random boost graph with no parallel edges. 
    /*!
     *  \param numVertices The number of vertices.
     *  \param numEdges The number of edges.
     *  \param minBlackWeight The lowest possible "black" edge weight.
     *  \param maxBlackWeight The maximum possible "black" edge weight.
     *  \param minRedWeight The lowest possible "red" edge weight.
     *  \param maxRedWeight The maximum possible "red" edge weight.
     *  
     *  - Two integer weights on each edge, called "black" and "red".
     *  - "black" weights are random integers chosen uniformly from 
     *    [minBlackWeight, maxBlackWeight].
     *  - "red" weights are random integers chosen uniformly from 
     *    [minRedWeight, maxRedWeight].
     *  
     *  \sa RandomGraphProblem and RandomGraphProblem()
     */
    void makeGraph(int numVertices, int numEdges, 
                   int minBlackWeight, int maxBlackWeight, 
                   int minRedWeight, int maxRedWeight);

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

    //! Return a reference to the underlying graph.
    Graph& graph();
    //! Return a reference to the source vertex (s).
    Vertex& source();
    //! Return a reference to the target vertex (t).
    Vertex& target();

  private:
    //! The underlying graph.
    Graph g_;
    //! The source vertex (s).
    Vertex s_;
    //! The target vertex (t).
    Vertex t_;
};


#endif  // EXAMPLE_CLASS_RANDOM_GRAPH_PROBLEM_H
