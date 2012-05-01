/*! \file SmallGraphProblem.h
 *  \brief Declaration of SmallGraphProblem, a simple 
 *         biobjective shortest path problem class used in BaseProblem.cpp
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef EXAMPLE_CLASS_SMALL_GRAPH_PROBLEM_H
#define EXAMPLE_CLASS_SMALL_GRAPH_PROBLEM_H


#include <string>
#include <vector>
#include <boost/graph/adjacency_list.hpp>

#include "../PointAndSolution.h"
#include "../BaseProblem.h"


using std::string;

using pareto_approximator::PointAndSolution;
using pareto_approximator::BaseProblem;


namespace small_graph_problem {


// Each edge has two weights: black and red.
class EdgeProperty 
{
  public:
    double black;
    double red;
};


typedef boost::adjacency_list<boost::listS, 
                              boost::vecS, 
                              boost::directedS, 
                              boost::no_property, 
                              EdgeProperty> Graph;
typedef Graph::vertex_descriptor Vertex;
typedef Graph::edge_descriptor   Edge;
typedef std::vector<Vertex>      PredecessorMap;


/*
 *  A small biobjective shortest path problem we'll use to test BaseProblem.
 *  
 *  Problem solutions will be s-t paths, where s is the problem's source 
 *  vertex (attribute s_) and t the problem's target vertex (attribute t_).
 *  
 *  Let P be an s-t path. We want solutions that minimize both of the 
 *  following objective functions (at the same time):
 *  - Black(P): the sum of all e.black weights, where e is an Edge on
 *    the path P.
 *  - Red(P): the sum of all e.red weights, where e is an Edge on the path P.
 */
class SmallGraphProblem : public BaseProblem<PredecessorMap>
{
  public:
    SmallGraphProblem();
    ~SmallGraphProblem();

    PointAndSolution<PredecessorMap> comb(std::vector<double>::const_iterator first, 
                                          std::vector<double>::const_iterator last);

    void makeGraph();


  private:
    // The problem's graph. (a boost adjacency list)
    Graph g_;
    // source vertex
    Vertex s_;
    // target vertex
    Vertex t_;
};


}  // namespace small_graph_problem


#endif  // EXAMPLE_CLASS_SMALL_GRAPH_PROBLEM_H
