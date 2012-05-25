/*! \file FloodVisitor.cpp
 *  \brief The implementation of the FloodVisitor class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <limits>

#include "biobjective_shortest_path_example_common.h"
#include "FloodVisitor.h"


using pareto_approximator::Point;
using pareto_approximator::NonDominatedSet;


namespace biobjective_shortest_path_example {


//! A simple constructor.
/*!
 *  \param source The source vertex of the graph.
 *  \param target The target vertex of the graph.
 *  \param numVertices The number of vertices in the graph.
 */
FloodVisitor::FloodVisitor(const Vertex & source, 
                           const Vertex & target, 
                           unsigned int numVertices) : source_(source), 
                                                       target_(target)
{
  vertexDistances_.assign(numVertices, NonDominatedSet<Point>());
}


void 
FloodVisitor::initializeVertex(const Vertex & u, const Graph &)
{
  if (u == source_)
    vertexDistances_[u].insert(Point(0.0, 0.0));
  else
    vertexDistances_[u].insert(Point(std::numeric_limits<double>::max(), 
                                              std::numeric_limits<double>::max()));
}


bool 
FloodVisitor::broadcastDistances(const Edge & e, const Graph & g)
{
  Vertex u = boost::source(e, g);
  Vertex v = boost::target(e, g);
  Point edgeWeight(g[e].black, g[e].red);

  bool insertedNewDistance = false;
  NonDominatedSet<Point>::iterator udi;
  for (udi = vertexDistances_[u].begin(); 
       udi != vertexDistances_[u].end(); ++udi) {
    insertedNewDistance |= vertexDistances_[v].insert(*udi + edgeWeight);
  }

  return insertedNewDistance;
}


NonDominatedSet<Point> 
FloodVisitor::getParetoPoints()
{
  return vertexDistances_[target_];
}


}  // namespace biobjective_shortest_path_example
