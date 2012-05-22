/*! \file main.cpp
 *  \brief The main program using RandomGraphProblem::computeConvexParetoSet() 
 *         to solve a biobjective shortest path problem.
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  We use a RandomGraphProblem object to represent a biobjective shortest 
 *  path problem on a random boost graph.
 *  
 *  All typedefs and class declarations are in RandomGraphProblem.h.
 *  
 *  \sa RandomGraphProblem and BaseProblem.
 */


#include <iostream>
#include <ctime>
#include <list>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include "../../PointAndSolution.h"
#include "RandomGraphProblem.h"


using std::cout;
using std::endl;


//! The example's main function.
/*!
 *  Will make a RandomGraphProblem instance and, if t is reachable, will 
 *  find a convex 0.001-approximate Pareto set and print its size, its 
 *  points and the problem solutions corresponding to those points.
 */
int 
main(void)
{
  using biobjective_shortest_path_example::RandomGraphProblem;
  using biobjective_shortest_path_example::PredecessorMap;

  // Initializations
  // =======================
  // Make a RandomGraphProblem instance with 1000 vertices and 100000 edges.
  // - "black" edge weights should be random integers in [1, 100] (uniformly)
  // - "red" edge weights should be random integers in [1, 100] (uniformly)
  // Reminder: All instances are created randomly so even instances with 
  // the same number of vertices and edges will probably be different.
  RandomGraphProblem rgp(1000, 100000, 1, 100, 1, 100);

  // Print problem info.
  cout << "Biobjective shortest path problem:" << endl
       << "- undirected random boost graph with no parallel edges" << endl
       << "- 1000 vertices and 100000 edges" << endl
       << "- two weights (\"black\" and \"red\") on each edge" << endl
       << "- \"black\" edge weights: random integers drawn uniformly from" 
       << " [1, 100]" << endl
       << "- \"red\" edge weights: random integers drawn uniformly from"
       << " [1, 100]" << endl
       << "- (random) source vertex s, (random) target vertex t" << endl
       << "- two objective functions to minimize: " << endl
       << "  (let P be an s-t path)" << endl
       << "  + Black(P): sum of \"black\" weights of all edges in P" << endl
       << "  + Red(P): sum of \"red\" weights of all edges in P" << endl
       << "- find a convex 1.001-approximation to the Pareto set" << endl
       << endl;

  // Check whether or not t is reachable (from s).
  if (rgp.isTargetReachable())
    cout << "Vertex t is reachable." << endl;
  else {
    // If it's not, don't bother trying to find shortest paths.
    cout << "Vertex t is not reachable! ";
    cout << "No point in continuing..." << endl << endl;
    return 1;
  }

  cout << "(computing... please wait a few seconds)" << endl << endl;
  // All the work (essentially 2 lines!)
  // =======================
  // Use RandomGraphProblem::computeConvexParetoSet() (inherited from 
  // BaseProblem) to find the approximate Pareto set.
  unsigned int numObjectives = 2;
  double approximationRatio = 0.001;
  std::list< PointAndSolution<PredecessorMap> > paretoSet;
  paretoSet = rgp.computeConvexParetoSet(numObjectives, approximationRatio);

  // Output
  // =======================
  cout << "A. (approximate) Pareto set size: " << paretoSet.size() << endl;
  cout << endl << "B. (approximate) Pareto (set) points: " << endl;
  std::list< PointAndSolution<PredecessorMap> >::iterator li;
  // Print each Pareto optimal point and the corresponding solution (path).
  for (li = paretoSet.begin(); li != paretoSet.end(); ++li) {
    cout << li->point << endl;
    rgp.printPath(li->solution);
  }

  return 0;
}


