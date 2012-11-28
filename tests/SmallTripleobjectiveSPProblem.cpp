/*! \file SmallTripleobjectiveSPProblem.cpp
 *  \brief Implementation of SmallTripleobjectiveSPProblem, a simple 
 *         triple-objective shortest path problem class used in 
 *         BaseProblemTest.cpp.
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <iostream>
#include <assert.h>
#include <map>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "SmallTripleobjectiveSPProblem.h"
#include "../Point.h"


using pareto_approximator::Point;


namespace small_tripleobjective_sp_problem {


SmallTripleobjectiveSPProblem::SmallTripleobjectiveSPProblem() 
{
  makeGraph();
}


SmallTripleobjectiveSPProblem::~SmallTripleobjectiveSPProblem() { }


void 
SmallTripleobjectiveSPProblem::makeGraph()
{
  s_ = boost::add_vertex(g_);
  Vertex m1 = boost::add_vertex(g_);
  Vertex m2 = boost::add_vertex(g_);
  Vertex m3 = boost::add_vertex(g_);
  Vertex m4 = boost::add_vertex(g_);
  t_ = boost::add_vertex(g_);

  bool ok;
  Edge sm1, sm2, sm3, sm4, m1t, m2t, m3t, m4t;
  boost::tie(sm1, ok) = boost::add_edge(s_, m1,  g_);
  boost::tie(sm2, ok) = boost::add_edge(s_, m2,  g_);
  boost::tie(sm3, ok) = boost::add_edge(s_, m3,  g_);
  boost::tie(sm4, ok) = boost::add_edge(s_, m4,  g_);
  boost::tie(m1t, ok) = boost::add_edge(m1, t_,  g_);
  boost::tie(m2t, ok) = boost::add_edge(m2, t_,  g_);
  boost::tie(m3t, ok) = boost::add_edge(m3, t_,  g_);
  boost::tie(m4t, ok) = boost::add_edge(m4, t_,  g_);
  
  g_[sm1].black = 1;
  g_[sm1].red   = 9;
  g_[sm1].blue  = 9;
  g_[m1t].black = 0;
  g_[m1t].red   = 0;
  g_[m1t].blue  = 0;
  g_[sm2].black = 9;
  g_[sm2].red   = 1;
  g_[sm2].blue  = 9;
  g_[m2t].black = 0;
  g_[m2t].red   = 0;
  g_[m2t].blue  = 0;
  g_[sm3].black = 9;
  g_[sm3].red   = 9;
  g_[sm3].blue  = 1;
  g_[m3t].black = 0;
  g_[m3t].red   = 0;
  g_[m3t].blue  = 0;
  g_[sm4].black = 5;
  g_[sm4].red   = 5;
  g_[sm4].blue  = 5;
  g_[m4t].black = 0;
  g_[m4t].red   = 0;
  g_[m4t].blue  = 0;
}


PointAndSolution<PredecessorMap> 
SmallTripleobjectiveSPProblem::comb(
                        std::vector<double>::const_iterator first, 
                        std::vector<double>::const_iterator last) const
{
  assert(last == first + 3);

  double blackWeight = *first;
  double redWeight = *(first + 1);
  double blueWeight = *(first + 2);

  std::map<Edge, double> weight;
  boost::graph_traits<Graph>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(g_); ei != ei_end; ++ei)
    weight[*ei] = blackWeight * g_[*ei].black + redWeight * g_[*ei].red + 
                  blueWeight * g_[*ei].blue;

  // weight map
  boost::associative_property_map< std::map<Edge, double> > w_map(weight);
  // predecessor map
  std::vector<Vertex> p_map(boost::num_vertices(g_));

  boost::dijkstra_shortest_paths(g_, s_, weight_map(w_map).predecessor_map(&p_map[0]));

  double blackDistance = 0;
  double redDistance = 0;
  double blueDistance = 0;
  Vertex v, w;
  w = t_;
  v = p_map[w];
  while (w != s_) {
    blackDistance += g_[boost::edge(v, w, g_).first].black;
    redDistance += g_[boost::edge(v, w, g_).first].red;
    blueDistance += g_[boost::edge(v, w, g_).first].blue;
    w = v;
    v = p_map[w];
  }

  return PointAndSolution<PredecessorMap>(Point(blackDistance, redDistance, 
                                                blueDistance), p_map);
}


}  // namespace small_tripleobjective_sp_problem
