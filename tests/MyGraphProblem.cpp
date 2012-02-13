#include <boost/graph/dag_shortest_paths.hpp>

#include "MyGraphProblem.h"
#include "../Point.h"


using pareto_approximator::Point;


MyGraphProblem::MyGraphProblem() 
{
  makeGraph();
}


void 
MyGraphProblem::makeGraph()
{
  s_ = boost::add_vertex(g_);
  Vertex u = boost::add_vertex(g_);
  Vertex d = boost::add_vertex(g_);
  Vertex l = boost::add_vertex(g_);
  Vertex m = boost::add_vertex(g_);
  Vertex r = boost::add_vertex(g_);
  t_ = boost::add_vertex(g_);

  boost::property_map<Graph, boost::vertex_name_t>::type vName = get(boost::vertex_name, g_);
  boost::put(vName, s_, 's');
  boost::put(vName, u, 'u');
  boost::put(vName, d, 'd');
  boost::put(vName, l, 'l');
  boost::put(vName, m, 'm');
  boost::put(vName, r, 'r');
  boost::put(vName, t_, 't');

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
  
  g_[su].name  = "su";
  g_[su].black = 7;
  g_[su].red   = 1;
  g_[sl].name  = "sl";
  g_[sl].black = 1;
  g_[sl].red   = 2;
  g_[sd].name  = "sd";
  g_[sd].black = 1;
  g_[sd].red   = 11;
  g_[um].name  = "um";
  g_[um].black = 1;
  g_[um].red   = 1;
  g_[ur].name  = "ur";
  g_[ur].black = 2;
  g_[ur].red   = 2;
  g_[ut].name  = "ut";
  g_[ut].black = 7;
  g_[ut].red   = 1;
  g_[dt].name  = "dt";
  g_[dt].black = 1;
  g_[dt].red   = 5;
  g_[lm].name  = "lm";
  g_[lm].black = 1;
  g_[lm].red   = 2;
  g_[ld].name  = "ld";
  g_[ld].black = 1;
  g_[ld].red   = 5;
  g_[mr].name  = "mr";
  g_[mr].black = 3;
  g_[mr].red   = 5;
  g_[md].name  = "md";
  g_[md].black = 1;
  g_[md].red   = 1;
  g_[rt].name  = "rt";
  g_[rt].black = 3;
  g_[rt].red   = 3;
  g_[rd].name  = "rd";
  g_[rd].black = 2;
  g_[rd].red   = 2;
}


PointAndSolution<PredecessorMap> 
MyGraphProblem::comb(double xWeight, double yWeight)
{
  boost::property_map<Graph, double EdgeProperty::*>::type   w_map   = boost::get(&EdgeProperty::weight,  g_);
  boost::property_map<Graph, boost::vertex_distance_t>::type d_map = boost::get(boost::vertex_distance, g_);
  std::vector<Vertex>                                        p_map(boost::num_vertices(g_));

  boost::graph_traits<Graph>::edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = edges(g_); ei != ei_end; ++ei)
    w_map[*ei] = xWeight * g_[*ei].black + yWeight * g_[*ei].red;

  boost::dag_shortest_paths(g_, s_, weight_map(w_map).distance_map(d_map).predecessor_map(&p_map[0]));

  double xDistance = 0;
  double yDistance = 0;
  Vertex v, w;
  w = t_;
  v = p_map[w];
  while (w != s_) {
    xDistance += g_[boost::edge(v, w, g_).first].black;
    yDistance += g_[boost::edge(v, w, g_).first].red;
    w = v;
    v = p_map[w];
  }

  return PointAndSolution<PredecessorMap>(Point(xDistance, yDistance), p_map);
}


