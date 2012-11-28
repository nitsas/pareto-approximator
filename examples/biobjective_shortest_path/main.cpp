/*! \file examples/biobjective_shortest_path/main.cpp
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
#include <cstdlib>
#include <time.h>
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



/*!
 *  \defgroup BiobjectiveShortestPathExample An example biobjective shortest path problem.
 *  
 *  @{
 */


//! The example's main function.
/*!
 *  Will make a RandomGraphProblem instance and, if t is reachable, will 
 *  find a convex Pareto set and print its size, its points and the 
 *  problem solutions corresponding to those points.
 */
int 
main(int argc, char * argv[])
{
  // Get the command line arguments (the random number generator's seed).
  int seed;
  if (argc > 2) {
    cout << "Too many arguments! Expected at most 1 integer argument "
         << "(a seed)." << endl;
    return(1);
  }
  else if (argc == 2) 
    // Use the input argument as a seed. 
    // Use 0 if the input argument is not an integer.
    seed = atoi(argv[1]);
  else 
    // Use the current time as a seed.
    seed = std::time(0);

  // Initializations
  // =========================================
  // Make a RandomGraphProblem instance with 1000 vertices and 100000 edges.
  // - "black" edge weights should be random integers in [1, 100] (uniformly)
  // - "red" edge weights should be random integers in [1, 100] (uniformly)
  // Reminder: All instances are created randomly so even instances with 
  // the same number of vertices and edges will almost surely be different.
  RandomGraphProblem rgp(1000, 100000, 1, 100, 1, 100, seed);

  // Print problem info.
  cout << "Biobjective shortest path problem:" << endl
       << "- undirected random boost graph with no parallel edges" << endl
       << "- random number generator's seed: " << seed << endl
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
       << "- find an eps-approximate convex Pareto set with eps = 1e-12" 
       << endl << endl;

  //cout << "Printing graph as a dot file (graph.dot)..." << endl;
  //rgp.printGraphToDotFile();
  //cout << "Done!" << endl << endl;

  // Check whether or not t is reachable (from s).
  if (rgp.isTargetReachable())
    cout << "Vertex t is reachable." << endl << endl;
  else {
    // If it's not, don't bother trying to find shortest paths.
    cout << "Vertex t is not reachable! ";
    cout << "No point in continuing..." << endl << endl;
    return 1;
  }

  cout << "(computing approximate convex Pareto set... "
       << "please wait a few seconds)" << endl << endl;
  // All the work (essentially 2 lines!)
  // =========================================
  // Use RandomGraphProblem::computeConvexParetoSet() (inherited from 
  // BaseProblem) to find the convex Pareto set.
  unsigned int numObjectives = 2;
  double approximationRatio = 1e-12;
  std::vector< PointAndSolution<PredecessorMap> > paretoSet;
  paretoSet = rgp.computeConvexParetoSet(numObjectives, approximationRatio);

  // Output (convex Pareto set)
  // =========================================
  cout << "A. approximate convex Pareto set size: " << paretoSet.size() << endl;
  cout << endl << "B. approximate convex Pareto set points: " << endl;
  cout << "   (they are all points of the exact Pareto set)" << endl;
  std::vector< PointAndSolution<PredecessorMap> >::iterator vi;
  // Print each Pareto optimal point and the corresponding solution (path).
  for (vi = paretoSet.begin(); vi != paretoSet.end(); ++vi) {
    cout << vi->point << endl;
    rgp.printPath(vi->solution);
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


/*! 
 *  @}
 */
