/* SmallGraphProblem.h */


#ifndef EXAMPLE_CLASS_SMALL_GRAPH_PROBLEM_H
#define EXAMPLE_CLASS_SMALL_GRAPH_PROBLEM_H


#include <string>
#include <boost/graph/adjacency_list.hpp>

#include "../PointAndSolution.h"
#include "../BaseProblem.h"


using std::string;

using pareto_approximator::PointAndSolution;
using pareto_approximator::BaseProblem;


class EdgeProperty 
{
  public:
    string name;
    double black;
    double red;
    double weight;
};


typedef boost::adjacency_list<boost::listS, 
                              boost::vecS, 
                              boost::directedS, 
                              boost::property<boost::vertex_name_t, char, 
                                boost::property<boost::vertex_distance_t, double> >,
                              EdgeProperty> Graph;
typedef Graph::vertex_descriptor Vertex;
typedef Graph::edge_descriptor   Edge;
typedef std::vector<Vertex>      PredecessorMap;


class SmallGraphProblem : public BaseProblem<PredecessorMap>
{
  public:
    SmallGraphProblem();
    ~SmallGraphProblem();

    PointAndSolution<PredecessorMap> comb(double xWeight, double yWeight);

    void makeGraph();


  private:
    Graph g_;
    Vertex s_;
    Vertex t_;
};


#endif  // EXAMPLE_CLASS_SMALL_GRAPH_PROBLEM_H
