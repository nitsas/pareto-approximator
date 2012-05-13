/*! \file SmallBiobjectiveSPProblem.cpp
 *  \brief Implementation of SmallBiobjectiveSPProblem, a simple 
 *         biobjective shortest path problem class used in BaseProblem.cpp
 *  \author Christos Nitsas
 *  \date 2012
 */


#include <iostream>
#include <assert.h>
#include <map>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/property_map/vector_property_map.hpp>

#include "SmallBiobjectiveSPProblem.h"
#include "../Point.h"


using pareto_approximator::Point;


namespace small_biobjective_sp_problem {


SmallBiobjectiveSPProblem::SmallBiobjectiveSPProblem() 
{
  makeGraph();
}


SmallBiobjectiveSPProblem::~SmallBiobjectiveSPProblem() { }


void 
SmallBiobjectiveSPProblem::makeGraph()
{
  s_ = boost::add_vertex(g_);
  Vertex u = boost::add_vertex(g_);
  Vertex d = boost::add_vertex(g_);
  Vertex l = boost::add_vertex(g_);
  Vertex m = boost::add_vertex(g_);
  Vertex r = boost::add_vertex(g_);
  t_ = boost::add_vertex(g_);

  bool ok;
  Edge su, sl, sd, um, ur, ut, dt, lm, ld, mr, md, rt, rd;
  boost::tie(su, ok) = boost::add_edge(s_, u,  g_);
  boost::tie(sl, ok) = boost::add_edge(s_, l,  g_);
  boost::tie(sd, ok) = boost::add_edge(s_, d,  g_);
  boost::tie(um, ok) = boost::add_edge(u,  m,  g_);
  boost::tie(ur, ok) = boost::add_edge(u,  r,  g_);
  boost::tie(ut, ok) = boost::add_edge(u,  t_, g_);
  boost::tie(dt, ok) = boost::add_edge(d,  t_, g_);
  boost::tie(lm, ok) = boost::add_edge(l,  m,  g_);
  boost::tie(ld, ok) = boost::add_edge(l,  d,  g_);
  boost::tie(mr, ok) = boost::add_edge(m,  r,  g_);
  boost::tie(md, ok) = boost::add_edge(m,  d,  g_);
  boost::tie(rt, ok) = boost::add_edge(r,  t_, g_);
  boost::tie(rd, ok) = boost::add_edge(r,  d,  g_);
  
  g_[su].black = 7;
  g_[su].red   = 1;
  g_[sl].black = 1;
  g_[sl].red   = 2;
  g_[sd].black = 1;
  g_[sd].red   = 11;
  g_[um].black = 1;
  g_[um].red   = 1;
  g_[ur].black = 2;
  g_[ur].red   = 2;
  g_[ut].black = 7;
  g_[ut].red   = 1;
  g_[dt].black = 1;
  g_[dt].red   = 5;
  g_[lm].black = 1;
  g_[lm].red   = 2;
  g_[ld].black = 1;
  g_[ld].red   = 5;
  g_[mr].black = 3;
  g_[mr].red   = 5;
  g_[md].black = 1;
  g_[md].red   = 1;
  g_[rt].black = 3;
  g_[rt].red   = 3;
  g_[rd].black = 2;
  g_[rd].red   = 2;
}


PointAndSolution<PredecessorMap> 
SmallBiobjectiveSPProblem::comb(std::vector<double>::const_iterator first, 
                                std::vector<double>::const_iterator last)
{
  assert(last == first + 2);

  double blackWeight = *first;
  double redWeight = *(first + 1);

  std::map<Edge, double> weight;
  boost::graph_traits<Graph>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(g_); ei != ei_end; ++ei)
    weight[*ei] = blackWeight * g_[*ei].black + redWeight * g_[*ei].red;

  // weight map
  boost::associative_property_map< std::map<Edge, double> > w_map(weight);
  // predecessor map
  std::vector<Vertex> p_map(boost::num_vertices(g_));
  // distance map
  boost::vector_property_map<double> d_map(boost::num_vertices(g_));

  // find all shortest paths starting from s_
  boost::bellman_ford_shortest_paths(g_, boost::num_vertices(g_), 
                                     weight_map(w_map).
                                     distance_map(d_map).
                                     predecessor_map(&p_map[0]).
                                     root_vertex(s_));

  double blackDistance = 0;
  double redDistance = 0;
  Vertex v, w;
  w = t_;
  v = p_map[w];
  while (w != s_) {
    blackDistance += g_[boost::edge(v, w, g_).first].black;
    redDistance += g_[boost::edge(v, w, g_).first].red;
    w = v;
    v = p_map[w];
  }

  return PointAndSolution<PredecessorMap>(Point(blackDistance, redDistance), p_map);
}


}  // namespace small_biobjective_sp_problem
