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
#include <list>
#include <boost/graph/adjacency_list.hpp>

#include "../../Point.h"
#include "../../PointAndSolution.h"
#include "../../NonDominatedSet.h"
#include "RandomGraphProblem.h"


using std::cout;
using std::endl;

using pareto_approximator::Point;
using pareto_approximator::PointAndSolution;
using pareto_approximator::NonDominatedSet;
using biobjective_shortest_path_example::RandomGraphProblem;
using biobjective_shortest_path_example::PredecessorMap;



//! The example's main function.
/*!
 *  Will make a RandomGraphProblem instance and, if t is reachable, will 
 *  find a convex 0.001-approximate Pareto set and print its size, its 
 *  points and the problem solutions corresponding to those points.
 */
int 
main(void)
{
  // Initializations
  // =========================================
  // Make a RandomGraphProblem instance with 1000 vertices and 100000 edges.
  // - "black" edge weights should be random integers in [1, 100] (uniformly)
  // - "red" edge weights should be random integers in [1, 100] (uniformly)
  // Reminder: All instances are created randomly so even instances with 
  // the same number of vertices and edges will almost surely be different.
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
       << "- find a convex Pareto set" << endl
       << endl;

  /*
  cout << "Printing graph as a dot file (graph.dot)..." << endl;
  rgp.printGraphToDotFile();
  cout << "Done!" << endl << endl;
  */

  // Check whether or not t is reachable (from s).
  if (rgp.isTargetReachable())
    cout << "Vertex t is reachable." << endl << endl;
  else {
    // If it's not, don't bother trying to find shortest paths.
    cout << "Vertex t is not reachable! ";
    cout << "No point in continuing..." << endl << endl;
    return 1;
  }

  cout << "(computing approximate convex Pareto set... please wait a few seconds)" << endl << endl;
  // All the work (essentially 2 lines!)
  // =========================================
  // Use RandomGraphProblem::computeConvexParetoSet() (inherited from 
  // BaseProblem) to find the approximate Pareto set.
  unsigned int numObjectives = 2;
  double approximationRatio = 0.0;
  std::list< PointAndSolution<PredecessorMap> > paretoSet;
  paretoSet = rgp.computeConvexParetoSet(numObjectives, approximationRatio);

  // Output (approximate convex Pareto set)
  // =========================================
  cout << "A. (approximate) convex Pareto set size: " << paretoSet.size() << endl;
  cout << endl << "B. (approximate) convex Pareto (set) points: " << endl;
  std::list< PointAndSolution<PredecessorMap> >::iterator li;
  // Print each Pareto optimal point and the corresponding solution (path).
  for (li = paretoSet.begin(); li != paretoSet.end(); ++li) {
    cout << li->point << endl;
    rgp.printPath(li->solution);
  }

  // Exact Pareto set
  // =========================================
  cout << endl << "(computing exact Pareto set... please wait a few seconds)" << endl << endl;
  NonDominatedSet<Point> exactParetoSet = rgp.findExactParetoSet();
  cout << "C. exact Pareto set size: " << exactParetoSet.size() << endl;
  cout << endl << "D. exact Pareto set points: " << endl;
  NonDominatedSet<Point>::iterator epsi;
  for (epsi = exactParetoSet.begin(); epsi != exactParetoSet.end(); ++epsi)
    cout << *epsi << endl;

  return 0;
}


