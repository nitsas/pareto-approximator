/*! \file RandomGraphProblem.cpp
 *  \brief The implementation of the RandomGraphProblem class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <map>
#include <list>
#include <time.h>
#include <boost/config.hpp>
#include <boost/graph/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/property_map/vector_property_map.hpp>

#include "RandomGraphProblem.h"
#include "../../Point.h"


using std::cout;
using std::endl;
using std::map;
using std::vector;

using boost::tie;
using boost::mt19937;
using boost::uniform_int;
using boost::variate_generator;
using boost::adjacency_list;
using boost::graph_traits;
using boost::num_vertices;

using pareto_approximator::Point;


//! A random number generator for uniformly distributed random integers.
/*!
 *  variate_generator template arguments:
 *  - mt19937 is a boost pseudorandom number generator (Mersenne twister)
 *  - uniform_int<> wraps mt19937 in a uniform integer distribution
 */
typedef variate_generator< mt19937&, uniform_int<> > UniformRandomIntGenerator;


//! Constructor. Make a biobjective shortest path problem instance.
/*!
 *  \param numVertices The number of vertices.
 *  \param numEdges The number of edges.
 *  \param minBlackWeight The lowest possible "black" edge weight. 
 *                        (must be \f$ \ge 0\f$ or the concepts of approximation
 *                        and of points dominating other points break down)
 *  \param maxBlackWeight The maximum possible "black" edge weight.
 *  \param minRedWeight The lowest possible "red" edge weight.
 *                      (must be \f$ \ge 0\f$ or the concepts of approximation
 *                      and of points dominating other points break down)
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
 *    (let P be an s-t path)
 *    + Black(P): the sum of "black" weights of all the edges in P.
 *    + Red(P): the sum of "red" weights of all the edges in P.
 *  
 *  \sa ~RandomGraphProblem() and makeGraph()
 */
RandomGraphProblem::RandomGraphProblem(int numVertices, int numEdges, 
                                       int minBlackWeight, int maxBlackWeight, 
                                       int minRedWeight, int maxRedWeight) 
{
  // Make a random graph.
  // - two weights on each edge
  // - "black" edge weights random integers in [minBlackWeight, maxBlackWeight]
  // - "red" edge weights random integers in [minRedWeight, maxRedWeight]
  makeGraph(numVertices, numEdges, minBlackWeight, maxBlackWeight, 
            minRedWeight, maxRedWeight);

  // Designate a source and a target vertex. 
  // - the first and last vertex created
  VertexIterator vi, vi_end;
  tie(vi, vi_end) = boost::vertices(g_);
  s_ = *vi;
  t_ = *(--vi_end);
}


//! Empty destructor.
/*!
 *  \sa RandomGraphProblem()
 */
RandomGraphProblem::~RandomGraphProblem() { }


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
void 
RandomGraphProblem::makeGraph(int numVertices, int numEdges, 
                              int minBlackWeight, int maxBlackWeight, 
                              int minRedWeight, int maxRedWeight)
{
  // Generate a random graph with numVertices vertices and numEdges edges.
  mt19937 generator(std::time(0));
  generate_random_graph(g_, numVertices, numEdges, generator);

  // The random integer weight distributions and generators.
  uniform_int<> blackDistribution(minBlackWeight, maxBlackWeight);
  uniform_int<> redDistribution(minRedWeight, maxRedWeight);
  UniformRandomIntGenerator randBlack(generator, blackDistribution);
  UniformRandomIntGenerator randRed(generator, redDistribution);

  //  Give each edge two weights ("black" and "red")
  //  - "black": an integer in [minBlackWeight, maxBlackWeight]
  //  - "red": an integer in [minRedWeight, maxRedWeight]
  EdgeIterator ei, ei_end;
  tie(ei, ei_end) = boost::edges(g_);
  while (ei != ei_end) {
    g_[*ei].black = randBlack();
    g_[*ei].red = randRed();
    ++ei;
  }
}


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
PointAndSolution<PredecessorMap> 
RandomGraphProblem::comb(double xWeight, double yWeight)
{
  map<Edge, double> weight;
  vector<Vertex>                                      p_map(num_vertices(g_));
  VertexIterator vi, vi_end;

  // Make the weight property map.
  EdgeIterator ei, ei_end;
  for (tie(ei, ei_end) = edges(g_); ei != ei_end; ++ei)
    weight[*ei] = xWeight * g_[*ei].black + yWeight * g_[*ei].red;
  boost::associative_property_map< map<Edge, double> > w_map(weight);

  // Find all shortest paths from s.
  boost::dijkstra_shortest_paths(g_, s_, weight_map(w_map).predecessor_map(&p_map[0]));

  double xDistance = 0;
  double yDistance = 0;
  Vertex v, w;
  w = t_;
  v = p_map[w];
  while (w != s_) {
    Edge e;
    bool ok;
    
    tie(e, ok) = boost::edge(v, w, g_);
    xDistance += g_[e].black;
    yDistance += g_[e].red;
    w = v;
    v = p_map[w];
  }

  return PointAndSolution<PredecessorMap>(Point(xDistance, yDistance), p_map);
}


//! Check if the target (t) is reachable.
/*!
 *  \return True iff there is at least one path that connects source (s) 
 *          and target (t).
 *  \return False iff there is no path that connects source (s) and 
 *          target (t).
 */
bool 
RandomGraphProblem::isTargetReachable()
{
  boost::vector_property_map<boost::default_color_type> c_map(num_vertices(g_));

  boost::breadth_first_search(g_, s_, color_map(c_map));
  if (c_map[t_] == boost::white_color) 
    return false;
  else
    return true;
}


//! A simple method that prints s-t paths.
/*!
 *  \param pred A map from each vertex to its predecessor in the path. 
 *              (Using a std::vector for the example.)
 *  
 *  Print the path as a series of vertex descriptors and edge weights 
 *  (formatted nicely).
 */
void 
RandomGraphProblem::printPath(const PredecessorMap& pred) const
{
  cout << "    ";
  if (pred[t_] == t_)
    // t was unreachable
    cout << "t was unreachable!" << endl;
  else {
    // there is a path from s to t, print it
    list< Edge > path;
    Vertex v;
    Edge e;
    bool exists;
  
    v = t_;
    while (v != s_) {
      tie(e, exists) = edge(pred[v], v, g_);
      path.push_front(e);
      v = pred[v];
    }

    list< Edge >::iterator li;
    cout << s_;
    for (li = path.begin(); li != path.end(); ++li)
      cout << "--(" << g_[*li].black << "," << g_[*li].red << ")-->" << boost::target(*li, g_);
    cout << endl;
  }
}


//! Return a reference to the underlying graph.
Graph& 
RandomGraphProblem::graph() 
{
  return g_;
}


//! Return a reference to the source vertex (s).
Vertex& 
RandomGraphProblem::source() 
{
  return s_;
}


//! Return a reference to the target vertex (t).
Vertex&
RandomGraphProblem::target() 
{
  return t_;
}


