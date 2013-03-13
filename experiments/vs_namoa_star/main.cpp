/*! \file experiments/vs_namoa_star/main.cpp
 *  \brief The main file running all the experiments vs the NAMOA\* algorithm.
 *  \author Christos Nitsas
 *  \date 2013
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <string>
#include <utility>
#include <cassert>
#include <limits>
#include <algorithm>
#include <sstream>
#include <sys/stat.h>

#include <Utilities/timer.h>

#include "experiments_vs_namoa_star_common.h"
#include "experiments_vs_namoa_star_utility.h"
#include "MultiobjectiveSpOnPmgProblem.h"


namespace pa = pareto_approximator;
namespace evns = experiments_vs_namoa_star;


/*! 
 *  \addtogroup ExperimentsVsNamoaStar Experiments vs the NAMOA Star algorithm.
 *  
 *  @{
 */


//! The experiments' main function.
int 
main(int argc, char * argv[])
{
  std::string mapPath = "/Users/chrisn/Programming/workspace/diplomatiki/dimacs-challenge-graphs/DIMACS9/";
  std::string mapName = "NY";
  std::string coordinatesFilename, graphFilename, distanceFilename, travelTimeFilename;
  std::string queriesFilePath = "/Users/chrisn/Programming/workspace/diplomatiki/pareto-approximator/experiments/vs_namoa_star/";
  std::string queriesFileName = "queries.txt";
  bool usingDimacs10Graph = false;

  // if there was a command line argument, assume it was the map name
  if (argc == 2) {
    mapName = std::string(argv[1]);
  }

  // is the map file a DIMACS-10 graph file?
  if ( (mapName == "belgium") || (mapName == "germany") || 
       (mapName == "italy") || (mapName == "luxembourg") || 
       (mapName == "netherlands") ) 
  {
    usingDimacs10Graph = true;
    mapPath = "/Users/chrisn/Programming/workspace/diplomatiki/dimacs-challenge-graphs/DIMACS10/";
  }

  if (usingDimacs10Graph) {
    // Using a DIMACS-10 graph.
    // make the filenames (from which we will read the graph)
    coordinatesFilename = mapPath + mapName + ".osm.xyz";
    graphFilename       = mapPath + mapName + ".osm.graph";

    // print the names of the graph files we will use
    std::cout << "Map files to be used: \n";
    std::cout << coordinatesFilename << "\n";
    std::cout << graphFilename << "\n" << std::endl;
  }
  else {
    // Using a DIMACS-9 graph.
    // make the filenames (from which we will read the graph)
    coordinatesFilename = mapPath + "USA-road-d." + mapName + ".co";
    distanceFilename    = mapPath + "USA-road-d." + mapName + ".gr";
    travelTimeFilename  = mapPath + "USA-road-t." + mapName + ".gr";

    // print the names of the graph files we will use
    std::cout << "Map files to be used: \n";
    std::cout << coordinatesFilename << "\n";
    std::cout << distanceFilename << "\n";
    std::cout << travelTimeFilename << "\n" << std::endl;
  }

  // print the name of the file containing the queries
  std::cout << "Queries file: \n" << (queriesFilePath + queriesFileName) << std::endl;

  // read (parse) the queries file
  std::list< std::pair<unsigned int, unsigned int> > queries;
  queries = evns::parseQueriesFile(queriesFilePath + queriesFileName);

  // print the queries we just read
  std::cout << "Read " << queries.size() 
            << (queries.size() == 1 ? " query:\n" : " queries: \n");
  std::list< std::pair<unsigned int, unsigned int> >::const_iterator qi;
  for (qi = queries.begin(); qi != queries.end(); ++qi) {
    std::cout << qi->first << " " << qi->second << "\n";
  }
  std::cout << std::endl;

  // make the problem instance (and read the graph)
  std::cout << "Reading graph:" << std::endl;
  evns::MultiobjectiveSpOnPmgProblem problem;
  if (usingDimacs10Graph) {
    problem.readDimacs10Graph(graphFilename, coordinatesFilename);
  }
  else {
    problem.readDimacs9Graph(distanceFilename, 
                             travelTimeFilename, 
                             coordinatesFilename);
  }

  // Make a directory where all the files containing points of Pareto 
  // and convex Pareto sets will be stored.
  mkdir("results", S_IRWXU | S_IRGRP | S_IROTH | S_IXGRP | S_IXOTH);

  unsigned int numObjectives = 2;        // CHANGE-HERE
  bool useNamoaStar;
  Timer timer;
  std::stringstream filename;
  std::cout << "\nRunning queries (" << mapName << " map, " << numObjectives << " objectives):\n";
  std::cout << "-----------------------------" << std::endl;
  for (qi = queries.begin(); qi != queries.end(); ++qi) {
    std::cout << "FROM NODE " << qi->first << " TO NODE " << qi->second << "\n";

    // Use the Chord algorithm with PGL's Dijkstra implementation.
    std::cout << "- Computing shortest path from node " << qi->first << " to node " 
              << qi->second << " using the Chord algorithm (& PGL's Dijkstra) ..." << std::endl;
    useNamoaStar = false;
    timer.start();
    std::vector< pa::PointAndSolution<evns::Path> > rcd = problem.runQuery(qi->first, qi->second, numObjectives, useNamoaStar);
    timer.stop();
    std::list< pa::PointAndSolution<evns::Path> > rcdlce = evns::computeLowerConvexEnvelopeOfPoints(rcd, numObjectives);
    rcd.assign(rcdlce.begin(), rcdlce.end());
    std::cout << "  Elapsed time: " << timer.getElapsedTime() << " sec, " << problem.getTimeSpentInComb() << " of them spent in comb(), " << problem.getNumCallsToComb() << " calls to comb()\n";
    std::cout << "  # Pareto points found: " << rcd.size() << "\n";
    // write points to file
    filename.str(std::string());
    filename << "results/query-" << qi->first << "-" << qi->second << "-Chord.txt";
    std::cout << "  Printing points to file " << filename.str() << " ..." << std::endl;
    evns::printPointsToFile(rcd.begin(), rcd.end(), filename.str());

    // clean node attributes
    problem.cleanNodeAttributes();

    // Use PGL's NAMOA*.
    std::cout << "- Computing shortest path from node " << qi->first << " to node "
              << qi->second << " using PGL's NAMOA* ..." << std::endl;
    useNamoaStar = true;
    timer.start();
    std::vector< pa::PointAndSolution<evns::Path> > rn = problem.runQuery(qi->first, qi->second, numObjectives, useNamoaStar);
    timer.stop();
    std::cout << "  Elapsed time: " << timer.getElapsedTime() << "\n";
    std::cout << "  # Pareto points found: " << rn.size() << "\n";
    // write points to file
    filename.str(std::string());
    filename << "results/query-" << qi->first << "-" << qi->second << "-NamoaStar.txt";
    std::cout << "  Printing points to file " << filename.str() << " ..." << std::endl;
    evns::printPointsToFile(rn.begin(), rn.end(), filename.str());

    // clean node attributes
    problem.cleanNodeAttributes();

    // Use the Chord algorithm with our A* implementation.
    std::cout << "- Computing shortest path from node " << qi->first << " to node " 
              << qi->second << " using the Chord algorithm (& our A*) ..." << std::endl;
    useNamoaStar = false;
    timer.start();
    std::vector< pa::PointAndSolution<evns::Path> > rca = problem.runQuery(qi->first, qi->second, numObjectives, useNamoaStar, true);
    timer.stop();
    std::list< pa::PointAndSolution<evns::Path> > rcalce = evns::computeLowerConvexEnvelopeOfPoints(rca, numObjectives);
    rca.assign(rcalce.begin(), rcalce.end());
    std::cout << "  Elapsed time: " << timer.getElapsedTime() << "\n";
    std::cout << "  # Pareto points found: " << rca.size() << "\n";
    
    // Compare Chord-with-our-A*'s results to Chord-with-PGL-Dijkstra's results.
    std::cout << "- Checking if Chord-with-our-A* found the same points as Chord-with-PGL-Dijkstra:\n";
    if ( (rca.size() == rcd.size()) && std::equal(rca.begin(), rca.end(), rcd.begin()) ) {
      std::cout << "  true\n";
    }
    else {
      std::cout << "  false\n";
      // write points to file
      filename.str(std::string());
      filename << "results/query-" << qi->first << "-" << qi->second << "-Chord-AStar.txt";
      std::cout << "  Printing points to file " << filename.str() << " ..." << std::endl;
      evns::printPointsToFile(rca.begin(), rca.end(), filename.str());
    }

    // clean node attributes
    problem.cleanNodeAttributes();

    // Compute the lower envelope of PGL's NAMOA*'s results and compare with Chord-with-PGL-Dijkstra's results.
    std::cout << "- Computing lower (convex) envelope of the set NAMOA* found ..." << std::endl;
    std::list< pa::PointAndSolution<evns::Path> > envelope = evns::computeLowerConvexEnvelopeOfPoints(rn, numObjectives);
    if (envelope.size() == rcd.size()) 
      std::cout << "  Has SAME number of points as Chord found." << "\n";
    else
      std::cout << "  Has " << envelope.size() << " points, i.e. " 
                << (envelope.size() > rcd.size() ? "MORE" : "LESS") << " than Chord found.\n";
    std::cout << "  Checking if the envelope contains all the points Chord found:\n";
    if ( std::includes(envelope.begin(), envelope.end(), rcd.begin(), rcd.end()) )
      std::cout << "  true\n";
    else
      std::cout << "  false\n";
    // write points to file
    filename.str(std::string());
    filename << "results/query-" << qi->first << "-" << qi->second << "-Envelope.txt";
    std::cout << "  Printing points to file " << filename.str() << " ..." << std::endl;
    evns::printPointsToFile(envelope.begin(), envelope.end(), filename.str());

    std::cout << "\n";
  }

  return 0;
}


/*! 
 *  @}
 */
