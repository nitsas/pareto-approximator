/*! \file FloodVisitor.h
 *  \brief The definition of the FloodVisitor class.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include "biobjective_shortest_path_example_common.h"
#include "../../Point.h"
#include "../../NonDominatedSet.h"


using pareto_approximator::Point;
using pareto_approximator::NonDominatedSet;


namespace biobjective_shortest_path_example {


//! A custom visitor for our flood algorithm.
/*!
 *  It floods the graph with 2-dimensional distances from the source s 
 *  and records all possible paths except paths that are definitely 
 *  not on the Pareto set (dominated paths).
 */
class FloodVisitor 
{
  public:
    //! A simple constructor.
    /*!
     *  \param source The source vertex of the graph.
     *  \param target The target vertex of the graph.
     *  \param numVertices The number of vertices in the graph.
     */
    FloodVisitor(const Vertex & source, const Vertex & target, 
                     unsigned int numVertices);

    NonDominatedSet<Point> getParetoPoints();

    void initializeVertex(const Vertex & u, const Graph & g);

    bool broadcastDistances(const Edge & e, const Graph & g);

  private:
    //! The source vertex. (for the shortest path problem)
    Vertex source_;

    //! The target vertex. (for the shortest path problem)
    Vertex target_;

    //! A vector of distances for each vertex. (distances from source_)
    /*!
     *  A vertex can have multiple distances, one for each path from 
     *  the source vertex. We will not be considering long paths 
     *  i.e. paths with cycles or paths that are definitely worse than 
     *  other paths we have already discovered.
     *
     *  We represent distances as Point instances.
     */
    std::vector< NonDominatedSet<Point> > vertexDistances_;
};


}  // namespace biobjective_shortest_path_example
