/*! \file SmallTripleobjectiveSPProblem.h
 *  \brief Declaration of SmallTripleobjectiveSPProblem, a simple 
 *         triple-objective shortest path problem class used in 
 *         BaseProblemTest.cpp.
 *  \author Christos Nitsas
 *  \date 2012
 */


#ifndef EXAMPLE_CLASS_SMALL_TRIPLEOBJECTIVE_SHORTEST_PATH_PROBLEM_H
#define EXAMPLE_CLASS_SMALL_TRIPLEOBJECTIVE_SHORTEST_PATH_PROBLEM_H


#include <string>
#include <vector>
#include <boost/graph/adjacency_list.hpp>

#include "../PointAndSolution.h"
#include "../BaseProblem.h"


using std::string;

using pareto_approximator::PointAndSolution;
using pareto_approximator::BaseProblem;


namespace small_tripleobjective_sp_problem {


// Each edge has three weights: black, red and blue.
class EdgeProperty 
{
  public:
    double black;
    double red;
    double blue;
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
 *  A small triple-objective shortest path problem we'll use to test 
 *  BaseProblem.
 *  
 *  Problem solutions will be s-t paths, where s is the problem's source 
 *  vertex (attribute s_) and t the problem's target vertex (attribute t_).
 *  
 *  Let P be an s-t path. We want solutions that minimize all the following
 *  objective functions (at the same time):
 *  - Black(P): the sum of e.black weights, over all edges e on the path P.
 *  - Red(P): the sum of e.red weights, over all edges e on the path P.
 *  - Blue(P): the sum of e.blue weights, over all edges e on the path P.
 */
class SmallTripleobjectiveSPProblem : public BaseProblem<PredecessorMap>
{
  public:
    SmallTripleobjectiveSPProblem();
    ~SmallTripleobjectiveSPProblem();

    PointAndSolution<PredecessorMap> comb(std::vector<double>::const_iterator first, 
                                          std::vector<double>::const_iterator last);

    void makeGraph();


  private:
    Graph g_;
    Vertex s_;
    Vertex t_;
};


}  // namespace small_tripleobjective_sp_problem


#endif  // EXAMPLE_CLASS_SMALL_TRIPLEOBJECTIVE_SHORTEST_PATH_PROBLEM_H
