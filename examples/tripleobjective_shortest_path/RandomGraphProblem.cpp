/*! \file examples/tripleobjective_shortest_path/RandomGraphProblem.cpp
 *  \brief The implementation of the RandomGraphProblem class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <assert.h>
#include <iostream>
#include <iterator>
#include <sstream>
#include <map>
#include <list>
#include <vector>
#include <queue>
#include <boost/config.hpp>
#include <boost/graph/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/property_map/vector_property_map.hpp>
#include <boost/graph/graphviz.hpp>

#include "RandomGraphProblem.h"
#include "FloodVisitor.h"


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
using pareto_approximator::NonDominatedSet;


/*!
 *  \addtogroup TripleobjectiveShortestPathExample An example tripleobjective shortest path problem.
 *  
 *  @{
 */


//! Everything needed for the example tripleobjective shortest path problem.
namespace tripleobjective_shortest_path_example {


//! A random number generator for uniformly distributed random integers.
/*!
 *  variate_generator template arguments:
 *  - mt19937 is a boost pseudorandom number generator (Mersenne twister)
 *  - uniform_int<> wraps mt19937 in a uniform integer distribution
 */
typedef variate_generator< mt19937&, uniform_int<> > UniformRandomIntGenerator;


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
RandomGraphProblem::RandomGraphProblem(int numVertices, int numEdges, 
                                       int minBlackWeight, int maxBlackWeight, 
                                       int minRedWeight, int maxRedWeight,
                                       int minGreenWeight, int maxGreenWeight, 
                                       int seed)
{
  // Make a random graph.
  // - three weights on each edge
  // - "black" edge weights random integers in [minBlackWeight, maxBlackWeight]
  // - "red" edge weights random integers in [minRedWeight, maxRedWeight]
  // - "green" edge weights random integers in [minGreenWeight, maxGreenWeight]
  makeGraph(numVertices, numEdges, minBlackWeight, maxBlackWeight, 
            minRedWeight, maxRedWeight, minGreenWeight, maxGreenWeight, seed);

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
void 
RandomGraphProblem::makeGraph(int numVertices, int numEdges, 
                              int minBlackWeight, int maxBlackWeight, 
                              int minRedWeight, int maxRedWeight, 
                              int minGreenWeight, int maxGreenWeight, 
                              int seed)
{
  // Generate a random graph with numVertices vertices and numEdges edges.
  mt19937 generator(seed);
  generate_random_graph(g_, numVertices, numEdges, generator);

  // The random integer weight distributions and generators.
  uniform_int<> blackDistribution(minBlackWeight, maxBlackWeight);
  uniform_int<> redDistribution(minRedWeight, maxRedWeight);
  uniform_int<> greenDistribution(minGreenWeight, maxGreenWeight);
  UniformRandomIntGenerator randBlack(generator, blackDistribution);
  UniformRandomIntGenerator randRed(generator, redDistribution);
  UniformRandomIntGenerator randGreen(generator, greenDistribution);

  //  Give each edge three weights ("black", "red" and "green")
  //  - "black": an integer in [minBlackWeight, maxBlackWeight]
  //  - "red": an integer in [minRedWeight, maxRedWeight]
  //  - "green": an integer in [minGreenWeight, maxGreenWeight]
  EdgeIterator ei, ei_end;
  tie(ei, ei_end) = boost::edges(g_);
  while (ei != ei_end) {
    g_[*ei].black = randBlack();
    g_[*ei].red = randRed();
    g_[*ei].green = randGreen();
    std::stringstream ss;
    ss << "(" << g_[*ei].black << ", " << g_[*ei].red << ", " 
       << g_[*ei].green << ")";
    g_[*ei].label = ss.str();
    ++ei;
  }
}


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
PointAndSolution<PredecessorMap> 
RandomGraphProblem::comb(std::vector<double>::const_iterator first, 
                         std::vector<double>::const_iterator last) const
{
  assert(std::distance(first, last) == 3);

  double xWeight, yWeight, zWeight;
  xWeight = *first;
  yWeight = *(first + 1);
  zWeight = *(first + 2);

  // weight property map
  map<Edge, double> weight;
  EdgeIterator ei, ei_end;
  for (tie(ei, ei_end) = edges(g_); ei != ei_end; ++ei)
    weight[*ei] = xWeight * g_[*ei].black + yWeight * g_[*ei].red + 
                  zWeight * g_[*ei].green;
  boost::associative_property_map< map<Edge, double> > w_map(weight);

  // predecessor property map
  vector<Vertex> p_map(num_vertices(g_));
  // distance property map
  boost::vector_property_map<double> d_map(boost::num_vertices(g_));

  // Find all shortest paths from s.
  boost::dijkstra_shortest_paths(g_, s_, weight_map(w_map).
                                         predecessor_map(&p_map[0]).
                                         distance_map(d_map));

  double xDistance = 0;
  double yDistance = 0;
  double zDistance = 0;
  Vertex v, w;
  w = t_;
  v = p_map[w];
  while (w != s_) {
    Edge e;
    bool ok;
    
    tie(e, ok) = boost::edge(v, w, g_);
    xDistance += g_[e].black;
    yDistance += g_[e].red;
    zDistance += g_[e].green;
    w = v;
    v = p_map[w];
  }

  Point point(xDistance, yDistance, zDistance);
  return PointAndSolution<PredecessorMap>(point, p_map);
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
  std::cout << "    ";
  if (pred[t_] == t_)
    // t was unreachable
    std::cout << "t was unreachable!" << std::endl;
  else {
    // there is a path from s to t, print it
    std::list< Edge > path;
    Vertex v;
    Edge e;
    bool exists;
  
    v = t_;
    while (v != s_) {
      tie(e, exists) = edge(pred[v], v, g_);
      path.push_front(e);
      v = pred[v];
    }

    std::list< Edge >::iterator li;
    std::cout << s_;
    for (li = path.begin(); li != path.end(); ++li)
      std::cout << "--(" << g_[*li].black << "," << g_[*li].red << ", "
                << g_[*li].green << ")-->" << boost::target(*li, g_);
    std::cout << std::endl;
  }
}


//! Find the exact Pareto set.
/*!
 *  We will just flood the graph with distances from the source starting
 *  from the source vertex.
 */
NonDominatedSet<Point> 
RandomGraphProblem::findExactParetoSet()
{
  using tripleobjective_shortest_path_example::FloodVisitor;

  // Flood the graph with distances from s_ (starting on s_).
  // We keep distances inside a FloodVisitor instance.
  std::queue<Vertex> q;
  FloodVisitor vis(s_, t_, boost::num_vertices(g_));
  
  VertexIterator vi, vi_end;
  for (tie(vi, vi_end) = boost::vertices(g_); vi != vi_end; ++vi) 
    vis.initializeVertex(*vi, g_);
  q.push(s_);
  while (q.size() > 0) {
    Vertex u = q.front(); q.pop();
    boost::graph_traits<Graph>::out_edge_iterator oei, oei_end;
    tie(oei, oei_end) = boost::out_edges(u, g_);
    for (; oei != oei_end; ++oei) {
      bool insertedAtLeastOne = vis.broadcastDistances(*oei, g_);
      Vertex v = boost::target(*oei, g_);
      assert(u != v);
      if (insertedAtLeastOne)
        q.push(v);
    }
  }

  return vis.getParetoPoints();
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


//! Print the graph to a dot (Graphviz) file.
void 
RandomGraphProblem::printGraphToDotFile(const char* filename)
{
  std::ofstream dotFile(filename);
  boost::dynamic_properties dp;
  dp.property("node_id", get(boost::vertex_index, g_));
  dp.property("label", get(&EdgeProperty::label, g_));
  write_graphviz_dp(dotFile, g_, dp);
}


}  // namespace tripleobjective_shortest_path_example


/*!
 *  @}
 */
