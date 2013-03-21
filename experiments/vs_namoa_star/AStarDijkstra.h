//! \file experiments/vs_namoa_star/AStarDijkstra.h
//! \brief The AStarDijkstra class. It is a simple implementation of the
//!        A\* algorithm wrapped inside a class.
//! \author Christos Nitsas
//! \date 2013
//!


#ifndef EXPERIMENTS_VS_NAMOA_STAR_A_STAR_DIJKSTRA_H
#define EXPERIMENTS_VS_NAMOA_STAR_A_STAR_DIJKSTRA_H


#include <assert.h>
#include <limits>
#include <cmath>

#include <Structs/Trees/priorityQueue.h>
//#include <Utilities/geographic.h>

#include "experiments_vs_namoa_star_common.h"
#include "experiments_vs_namoa_star_utility.h"


//! 
//! \addtogroup ExperimentsVsNamoaStar Experiments vs the NAMOA Star algorithm.
//! 
//! @{
//!


//! Everything needed for the experiments vs the NAMOA\* algorithm.
namespace experiments_vs_namoa_star {


/*
//! \brief Always returns 0 for every node.
//!
//! Use of this heuristic makes A* equivalent to the simple Dijkstra 
//! algorithm.
//!
template <class GraphType> 
class BlindHeuristic
{
  public:
    //! \brief Iterator to the underlying graph's nodes.
    typedef typename GraphType::NodeIterator NodeIterator;

    //! \brief Iterator to the underlying graph's edges.
    typedef typename GraphType::EdgeIterator EdgeIterator;

    //! \brief Constructor. (empty)
    //!
    //! \param graph The graph whose nodes we will be asked to estimate 
    //!        distances for.
    //!
    BlindHeuristic(GraphType & graph) : graph_(graph)
    {
    }

    //! \brief Compute an estimate of the distance between the given nodes. 
    //!
    //! Always returns zero.
    //!
    double 
    getHeuristicValue(const NodeIterator &, const NodeIterator &) const 
    {
      return 0;
    }

  private:
    //! The graph whose nodes we will be asked to estimate distances for.
    GraphType & graph_;
};
*/


//! \brief A simple heuristic using the Euclidean distance between nodes.
//!        (for DIMACS 10 graphs)
//!
//! It uses the nodes' x and y attributes (coordinates).
//!
//! It sets the nodes' "heuristicList" attributes (one value for every 
//! objective).
//!
template <class GraphType>
class EuclideanHeuristic
{
  public:
    //! \brief Iterator to the underlying graph's nodes.
    typedef typename GraphType::NodeIterator NodeIterator;

    //! \brief Iterator to the underlying graph's edges.
    typedef typename GraphType::EdgeIterator EdgeIterator;

    //! \brief Constructor. (empty)
    //!
    //! \param graph The graph whose nodes we will be asked to estimate 
    //!        distances for.
    //!
    EuclideanHeuristic(GraphType & graph) : graph_(graph), 
                                            maxSpeed_(0.0)
//                                            , maxEdgeDistance_(0) // CHANGE-HERE
    {
      NodeIterator u, v, lastNode;
      EdgeIterator e, lastEdge;
      double speed;

      for (u = graph_.beginNodes(), lastNode = graph_.endNodes(); 
           u != lastNode; ++u)
      {
        for (e = graph_.beginEdges(u), lastEdge = graph_.endEdges(u); 
             e != lastEdge; ++e)
        {
          v = graph_.target(e);

          if (e->criteriaList[0] == 0 && e->criteriaList[1] == 0) {
            continue;
          }

          // temporary?
          assert(e->criteriaList[1] != 0);

          speed = double(e->criteriaList[0]) / e->criteriaList[1];
          if (speed > maxSpeed_)
            maxSpeed_ = speed;

          // CHANGE-HERE
          /*
          if (e->criteriaList[0] > maxEdgeDistance_) {
            maxEdgeDistance_ = e->criteriaList[0];
          }
          */
        }
      }
    }

    //! \brief Initialize every node's Node::heuristicList attribute.
    //!        (one heuristic value for every objective)
    //!
    //! \param target The target node. Every node u's heuristic values are 
    //!        estimates for the cost of the shortest path between u and t.
    //!
    void 
    initHeuristicLists(const NodeIterator & target)
    {
      NodeIterator u, lastNode;

      for (u = graph_.beginNodes(), lastNode = graph_.endNodes(); 
           u != lastNode; ++u)
      {
        // first, the heuristic for the "Distance" objective:
        u->heuristicList[0] = euclideanDistance(u->x, u->y, target->x, target->y);
        // next, the heuristic for the "Travel Time" objective:
        u->heuristicList[1] = u->heuristicList[0] / maxSpeed_;
        // CHANGE-HERE
        /*
        // next, the heuristic for the "Number of Hops" objective:
        u->heuristicList[2] = u->heuristicList[0] / maxEdgeDistance_;
        */
      }
    }

  private:
    //! The graph whose nodes we will be asked to estimate distances for.
    GraphType & graph_;

    //! The max speed over any edge on the graph.
    double maxSpeed_;

    // CHANGE-HERE
    /*
    //! The max edge distance.
    double maxEdgeDistance_;
    */
};


//! \brief A simple heuristic depending on the Great Circle distance 
//!        between nodes.
//!
//! The Great Circle distance takes into account that the Earth is not flat.
//! The Great Circle distance is the shortest distance between any two 
//! points on the surface of a sphere, measured along a path on the 
//! surface of the sphere (as opposed to going through the sphere's 
//! interior).
//! 
//! It uses the nodes' x and y attributes (coordinates).
//! 
//! Obviously, we assume that the nodes' x and y attributes represent 
//! the nodes' latitude and longitude.
//!
//! This uses the greatCircleDistance() function, which computes Great 
//! Circle distances in meters.
//!
template <class GraphType>
class GreatCircleDistanceHeuristic
{
  public:
    //! \brief Iterator to the underlying graph's nodes.
    typedef typename GraphType::NodeIterator NodeIterator;

    //! \brief Iterator to the underlying graph's edges.
    typedef typename GraphType::EdgeIterator EdgeIterator;

    //! \brief Constructor. (empty)
    //!
    //! \param graph The graph whose nodes we will be asked to estimate 
    //!        distances for.
    //!
    GreatCircleDistanceHeuristic(GraphType & graph) : graph_(graph), 
                                                      maxSpeed_(0.0)
//                                                      , maxEdgeDistance_(0) // CHANGE-HERE
    {
      NodeIterator u, v, lastNode;
      EdgeIterator e, lastEdge;
      double speed;

      for (u = graph_.beginNodes(), lastNode = graph_.endNodes(); 
           u != lastNode; ++u)
      {
        for (e = graph_.beginEdges(u), lastEdge = graph_.endEdges(u); 
             e != lastEdge; ++e)
        {
          v = graph_.target(e);

          if (e->criteriaList[0] == 0 && e->criteriaList[1] == 0) {
            continue;
          }

          // temporary?
          assert(e->criteriaList[1] != 0);

          speed = double(e->criteriaList[0]) / e->criteriaList[1];
          if (speed > maxSpeed_)
            maxSpeed_ = speed;

          // CHANGE-HERE
          /*
          if (e->criteriaList[0] > maxEdgeDistance_) {
            maxEdgeDistance_ = e->criteriaList[0];
          }
          */
        }
      }
    }

    //! \brief Initialize every node's Node::heuristicList attribute.
    //!        (one heuristic value for every objective)
    //!
    //! \param target The target node. Every node u's heuristic values are 
    //!        estimates for the cost of the shortest path between u and t.
    //!
    void 
    initHeuristicLists(const NodeIterator & target)
    {
      NodeIterator u, lastNode;

      for (u = graph_.beginNodes(), lastNode = graph_.endNodes(); 
           u != lastNode; ++u)
      {
        // first, the heuristic for the "Distance" objective:
        u->heuristicList[0] = floor(greatCircleUnderestimate(double(u->x)/100000, double(u->y)/100000, 
                                                       double(target->x)/100000, double(target->y)/100000));
        // next, the heuristic for the "Travel Time" objective:
        u->heuristicList[1] = u->heuristicList[0] / maxSpeed_;
        // CHANGE-HERE
        /*
        // next, the heuristic for the "Number of Hops" objective:
        u->heuristicList[2] = u->heuristicList[0] / maxEdgeDistance_;
        */
      }
    }

  private:
    //! The graph whose nodes we will be asked to estimate distances for.
    GraphType & graph_;

    //! The max speed over any edge on the graph.
    double maxSpeed_;

    // CHANGE-HERE
    /*
    //! The max edge distance.
    unsigned int maxEdgeDistance_;
    */
};


//! \brief A class containing a simple implementation of the A\* algorithm.
//!
template <class GraphType> 
class AStarDijkstra
{
  public:
    //! \brief Iterator to the underlying graph's nodes.
    typedef typename GraphType::NodeIterator         NodeIterator;

    //! \brief Node descriptor for the underlying graph's nodes.
    typedef typename GraphType::NodeDescriptor       NodeDescriptor;

    //! \brief Iterator to the underlying graph's edges.
    typedef typename GraphType::EdgeIterator         EdgeIterator;

    //! \brief Edge weight type.
    typedef double                                   WeightType;

    //! \brief The type of the priority we will use for the algorithm.
    typedef PriorityQueue<WeightType, 
                          NodeIterator, HeapStorage> PriorityQueueType;

    //! \brief Simple constructor.
    //! 
    //! \param graph The graph to run the algorithm on.
    //! \param timestamp A pointer to a global timestamp or at least to 
    //!        a variable containing a timestamp greater than or equal to 
    //!        the maximum of the graph's nodes' timestamps. The timestamp 
    //!        does not have to be an actual system timestamp, just an 
    //!        unsigned integer. 
    //!
    //! The given timestamp should be greater than or equal to the maximum 
    //! of the nodes' timestamps (Node::timestamp attribute). The 
    //! algorithm's timestamp represents the last "time" a query was run 
    //! (or the "time" the algorithm object was initialized, if no query 
    //! has been run yet).
    //! 
    //! Reminder: Each node's timestamp represents the last "time" an 
    //!           operation was performed on the node. 
    //! 
    AStarDijkstra(GraphType & graph, unsigned int * timestamp) 
                 : graph_(graph), timestamp_(timestamp)
    {
    }

    //! \brief Run a source-target query and return the path's distance.
    //!
    //! \param source (An iterator to) The source node.
    //! \param target (An iterator to) The target node.
    //! \return The cost/weight of the path from the source to the target.
    //! 
    //! The actual path can be reconstructed, if need be, by following the 
    //! nodes' predecessor pointers (Node::pred attributes) starting at the 
    //! target node.
    //!
    WeightType runQuery(const NodeIterator & source, 
                        const NodeIterator & target)
    {
      NodeIterator u, v, lastNode;
      EdgeIterator e, lastEdge;
      WeightType tempFScore;
      PriorityQueueType openNodesQueue;

      // increment the timestamp; now every node's timestamp should be 
      // lower than the algorithm's
      // Reminder: We consider nodes with timestamps lower than the 
      //           algorithm's timestamp uninitialized.
      ++(*timestamp_);
      // initialize the source node's (relevant) attributes
      source->dist = 0.0;
      source->timestamp = *timestamp_;
      source->pred = graph_.nilNodeDescriptor();

      // compute the source node's fScore and insert it in the OPEN "list"
      source->fScore = source->dist + source->heuristicValue;
      openNodesQueue.insert(source->fScore, source, &(source->pqitem));
      source->closed = false;

      // while there are nodes in the OPEN "list"
      while (not openNodesQueue.empty()) {
        // remove the node with minimum fScore (i.e. distance from 
        // source node + heuristic estimate to target node) from the 
        // OPEN "list" and put it in the CLOSED "list"
        u = openNodesQueue.minItem();
        openNodesQueue.popMin();
        u->closed = true;

        if (u == target) 
          // since the heuristic is admissible (we assumed that), no other 
          // path will be shorter than the one that got us to t now
          break;

        // for each successor of u:
        for (e = graph_.beginEdges(u), lastEdge = graph_.endEdges(u); 
             e != lastEdge; ++e) {
          v = graph_.target(e);

          // compute the node's fScore if we were to pass through u
          tempFScore = u->dist + e->weight + v->heuristicValue;

          assert(v->timestamp <= *timestamp_);
          if (v->timestamp < *timestamp_) {
            // we consider the node uninitialized, we ignore its "closed"
            // attribute (it is not valid); we know that the node is in 
            // neither the OPEN or the CLOSED "list"
            
            // initialize the node's (relevant) attributes
            v->pred = u->getDescriptor();
            v->dist = u->dist + e->weight;
            v->fScore = tempFScore;
            v->timestamp = *timestamp_;

            // insert the node in the OPEN "list"
            openNodesQueue.insert(v->fScore, v, &(v->pqitem));
            v->closed = false;
          }
          else if (openNodesQueue.contains(&v->pqitem) && tempFScore < v->fScore) {
            // we have encountered the node before (it has been initialized) 
            // so we checked if it is in the OPEN "list" (openNodesQueue); 
            // the node is already in the OPEN "list"; we have found a 
            // better path to the node

            // update the node's (relevant) attributes
            v->pred = u->getDescriptor();
            v->dist = u->dist + e->weight;
            v->fScore = tempFScore;

            // update (decrease) the node's fScore in the OPEN "list"
            openNodesQueue.decrease(v->fScore, &(v->pqitem));
          }
          else if (v->closed && tempFScore < v->fScore) {
            // we have encountered the node before (it has been initialized) 
            // so its "closed" attribute is valid; the node is in the 
            // CLOSED "list"; we have found a better path to the node

            // update the node's (relevant) attributes
            v->pred = u->getDescriptor();
            v->dist = u->dist + e->weight;
            v->fScore = tempFScore;

            // reinsert the node in the OPEN "list" (it is removed from CLOSED)
            openNodesQueue.insert(v->fScore, v, &(v->pqitem));
            v->closed = false;
          }
          // else we have already found a better path to the node; 
          // we ignore the current path

        }  // end for
      }  // end while

      // target's dist attribute has already been initialized 
      // (when it was inserted in the OPEN "list")

      return target->dist;
    }

    //! \brief Set every node's relevant attributes to default/initial values.
    //!
    //! Relevant attributes are Node::timestamp, Node::pred, Node::dist, 
    //! Node::fScore, Node::closed and Node::pqitem.
    //!
    //! The algorithm also uses the nodes' Node::x and Node::y attributes 
    //! but they are the node's coordinates, we won't change them.
    //!
    void cleanRelevantNodeAttributes()
    {
      NodeIterator ni, lastNode;
      for (ni = graph_.beginNodes(), lastNode = graph_.endNodes(); 
           ni != lastNode; ++ni)
      {
        // set the node's last used time to be the same as the 
        // algorithm's last run query time
        ni->timestamp = *timestamp_;
        // nil predecessor
        ni->pred = graph_.nilNodeDescriptor();
        // maximum distance
        ni->dist = std::numeric_limits<double>::infinity();
        // maximum fScore
        ni->fScore = std::numeric_limits<double>::infinity();
        // not in the CLOSED "list"
        ni->closed = false;
        // invalid value for the "pointer" to node's position in the OPEN "list"
        ni->pqitem = std::numeric_limits<unsigned int>::max();
      }
    }

  private:
    //! The graph the algorithm will run on.
    GraphType & graph_;

    //! \brief A timestamp. 
    //!
    //! The algorithm will increment it at the start of each query. Nodes 
    //! with timestamps lower than the algorithm's timestamp are considered 
    //! uninitialized.
    //!
    unsigned int * timestamp_;
};


}  // namespace experiments_vs_namoa_star


//!
//! @}
//!


#endif  // EXPERIMENTS_VS_NAMOA_STAR_A_STAR_DIJKSTRA_H
