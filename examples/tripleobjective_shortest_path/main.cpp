/*! \file main.cpp
 *  \brief The main program using RandomGraphProblem::computeConvexParetoSet() 
 *         to solve a tripleobjective shortest path problem.
 *  \author Christos Nitsas
 *  \date 2012
 *  
 *  We use a RandomGraphProblem object to represent a tripleobjective 
 *  shortest path problem on a random boost graph.
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
using tripleobjective_shortest_path_example::RandomGraphProblem;
using tripleobjective_shortest_path_example::PredecessorMap;



/*!
 *  \defgroup TripleobjectiveShortestPathExample An example tripleobjective shortest path problem.
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
  // Make a RandomGraphProblem instance with 100 vertices and 800 edges.
  // - "black" edge weights should be random integers in [1, 100] (uniformly)
  // - "red" edge weights should be random integers in [1, 100] (uniformly)
  // - "green" edge weights should be random integers in [1, 100] (uniformly)
  // Reminder: All instances are created randomly so even instances with 
  // the same number of vertices and edges will almost surely be different.
  RandomGraphProblem rgp(100, 800, 1, 100, 1, 100, 1, 100, seed);

  // Print problem info.
  cout << "Triple-objective shortest path problem:" << endl
       << "- undirected random boost graph with no parallel edges" << endl
       << "- random number generator's seed: " << seed << endl
       << "- 100 vertices and 800 edges" << endl
       << "- three weights (\"black\", \"red\" and \"green\") on each edge" 
       << endl
       << "- \"black\" edge weights: random integers drawn uniformly from" 
       << " [1, 100]" << endl
       << "- \"red\" edge weights: random integers drawn uniformly from"
       << " [1, 100]" << endl
       << "- \"green\" edge weights: random integers drawn uniformly from"
       << " [1, 100]" << endl
       << "- three objective functions to minimize: " << endl
       << "  (let P be an s-t path)" << endl
       << "  + Black(P): sum of \"black\" weights of all edges in P" << endl
       << "  + Red(P): sum of \"red\" weights of all edges in P" << endl
       << "  + Green(P): sum of \"green\" weights of all edges in P" << endl
       << "- find a convex Pareto set" << endl
       << endl;

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

  cout << "(computing convex Pareto set... please wait a few seconds)" << endl << endl;
  // All the work (essentially 2 lines!)
  // =========================================
  // Use RandomGraphProblem::computeConvexParetoSet() (inherited from 
  // BaseProblem) to find the Pareto set.
  unsigned int numObjectives = 3;
  double approximationRatio = 1e-12;
  std::list< PointAndSolution<PredecessorMap> > paretoSet;
  paretoSet = rgp.computeConvexParetoSet(numObjectives, approximationRatio);

  // Output (convex Pareto set)
  // =========================================
  cout << "A. convex Pareto set size: " << paretoSet.size() << endl;
  cout << endl << "B. convex Pareto (set) points: " << endl;
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


/*! 
 *  @}
 */


