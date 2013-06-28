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
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <csignal>
#include <cstdio>
#include <cerrno>

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


//! \brief A handler for the SIGALRM signal. (empty)
//!
void 
sigalrmHandler(int sig)
{
  // empty, we just need the alarm to un-block the waitpid() function
}


//! \brief Fork and run a query.
//!
void 
forkAndRunQuery(evns::MultiobjectiveSpOnPmgProblem & problem, 
                unsigned int sourceNodeId, unsigned int targetNodeId, 
                struct itimerval timeout, unsigned int numObjectives, 
                bool useNamoaStar, bool useAStar=false)
{
  std::stringstream filename;
  Timer timer;

  // now fork
  pid_t pid = fork();
  if (pid < 0) {
    // could not fork()
    std::cerr << "Failed to fork.. Exiting" << std::endl;
    exit(1);
  }
  else if (pid == 0) {
    // child process

    if (useNamoaStar) {
      // Flags: useNamoaStar is true, useAStar is irrelevant
      // Use PGL's NAMOA*.
      std::cout << "- Computing shortest path from node " << sourceNodeId << " to node "
                << targetNodeId << " using PGL's NAMOA* ..." << std::endl;
      timer.start();
      std::vector< pa::PointAndSolution<evns::Path> > rn = problem.runQuery(sourceNodeId, targetNodeId, numObjectives, useNamoaStar);
      timer.stop();
      std::cout << "  Elapsed time: " << timer.getElapsedTime() << "\n";
      std::cout << "  # Pareto points found: " << rn.size() << "\n";
      // write points to file
      filename.str(std::string());
      filename << "results/query-" << sourceNodeId << "-" << targetNodeId << "-NamoaStar.txt";
      std::cout << "  Printing points to file " << filename.str() << " ..." << std::endl;
      evns::printPointsToFile(rn.begin(), rn.end(), filename.str());

      // Compute the convex hull of PGL's NAMOA*'s results and compare with Chord-with-PGL-Dijkstra's results.
      std::list< pa::PointAndSolution<evns::Path> > hull = pa::utility::computeConvexHull(rn, numObjectives);
      std::cout << "- Comparing the convex hull of the set NAMOA* found with the set Chord found ..." << std::endl;
      // read the points Chord found (from the output file)
      filename.str(std::string());
      filename << "results/query-" << sourceNodeId << "-" << targetNodeId << "-Chord.txt";
      std::vector< pa::PointAndSolution<evns::Path> > rcd = evns::readPointsFromFile(filename.str());
      // CHANGE--HERE
      if (hull.size() == rcd.size()) 
        std::cout << "  Has SAME number of points as Chord found." << "\n";
      else {
        std::cout << "  Has " << hull.size() << " points, i.e. " << (hull.size() > rcd.size() ? "MORE" : "LESS") 
                  << " than Chord found. (Chord found " << rcd.size() << ")\n";
      }
      std::cout << "  Checking if the convex hull contains all the points Chord found:\n";
      if ( std::includes(hull.begin(), hull.end(), rcd.begin(), rcd.end()) )
        std::cout << "  true\n";
      else {
        std::cout << "  false";
        if ( std::includes(rn.begin(), rn.end(), rcd.begin(), rcd.end()) )
          std::cout << " (the complete Pareto set does though)\n";
        else
          std::cout << " (neither does the complete Pareto set)\n";
      }
      // write the convex hull's points to the output
      filename.str(std::string());
      filename << "results/query-" << sourceNodeId << "-" << targetNodeId << "-ConvexHull.txt";
      std::cout << "- Computing the convex hull of the set NAMOA* found ...\n"
                << "  Printing points to file " << filename.str() << " ..." << std::endl;
      evns::printPointsToFile(hull.begin(), hull.end(), filename.str());
      // Compare the convex hull of PGL's NAMOA*'s results with Chord-with-A*'s results.
      std::cout << "- Comparing the convex hull of the set NAMOA* found with the set Chord-with-A* found ..." << std::endl;
      // read the points Chord-with-A* found (from the output file)
      filename.str(std::string());
      filename << "results/query-" << sourceNodeId << "-" << targetNodeId << "-Chord-AStar.txt";
      std::vector< pa::PointAndSolution<evns::Path> > rca = evns::readPointsFromFile(filename.str());
      // CHANGE--HERE
      if (hull.size() == rca.size()) 
        std::cout << "  Has SAME number of points as Chord-with-A* found.\n";
      else {
        std::cout << "  Has " << hull.size() << " points, i.e. " << (hull.size() > rca.size() ? "MORE" : "LESS") 
                  << " than Chord-with-A* found. (Chord-with-A* found " << rca.size() << ")\n";
      }
      std::cout << "  Checking if the convex hull contains all the points Chord-with-A* found:\n";
      if ( std::includes(hull.begin(), hull.end(), rca.begin(), rca.end()) )
        std::cout << "  true" << std::endl;
      else {
        std::cout << "  false";
        if ( std::includes(rn.begin(), rn.end(), rca.begin(), rca.end()) )
          std::cout << " (the complete Pareto set does though)" << std::endl;
        else
          std::cout << " (neither does the complete Pareto set)" << std::endl;
      }
    }
    else if (not useAStar) {
      // Flags: useNamoaStar is false, useAStar is false
      // Use Chord with PGL's Dijkstra implementation
      std::cout << "- Computing shortest path from node " << sourceNodeId << " to node " 
                << targetNodeId << " using the Chord algorithm (& PGL's Dijkstra) ..." << std::endl;
      timer.start();
      std::vector< pa::PointAndSolution<evns::Path> > rcd = problem.runQuery(sourceNodeId, targetNodeId, numObjectives, useNamoaStar);
      timer.stop();
      std::cout << "  (Shortest) distance between nodes: " << rcd[0].point[0] << "\n";
      std::cout << "  Elapsed time: " << timer.getElapsedTime() << " sec, " << problem.getTimeSpentInComb() << " of them spent in comb(), " << problem.getNumCallsToComb() << " calls to comb()\n";
      std::list< pa::PointAndSolution<evns::Path> > rcdch = pa::utility::computeConvexHull(rcd, numObjectives);
      std::cout << "  # Pareto points found: " << rcdch.size() << "\n";
      // write points to file
      filename.str(std::string());
      filename << "results/query-" << sourceNodeId << "-" << targetNodeId << "-Chord.txt";
      std::cout << "  Printing points to file " << filename.str() << " ..." << std::endl;
      evns::printPointsToFile(rcdch.begin(), rcdch.end(), filename.str());
    }
    else {
      // Flags: useNamoaStar is false, useAStar is true
      // Use Chord with PGL's A* implementation
      std::cout << "- Computing shortest path from node " << sourceNodeId << " to node " 
                << targetNodeId << " using the Chord algorithm (& PGL's A*) ..." << std::endl;
      timer.start();
      std::vector< pa::PointAndSolution<evns::Path> > rca = problem.runQuery(sourceNodeId, targetNodeId, numObjectives, useNamoaStar, useAStar);
      timer.stop();
      std::cout << "  (Shortest) distance between nodes: " << rca[0].point[0] << "\n";
      std::cout << "  Elapsed time: " << timer.getElapsedTime() << " sec, " << problem.getTimeSpentInComb() << " of them spent in comb(), " << problem.getNumCallsToComb() << " calls to comb()\n";
      std::list< pa::PointAndSolution<evns::Path> > rcach = pa::utility::computeConvexHull(rca, numObjectives);
      std::cout << "  # Pareto points found: " << rcach.size() << "\n";
      
      // Compare Chord-with-our-A*'s results to Chord-with-PGL-Dijkstra's results.
      std::cout << "- Checking if Chord-with-PGL-A* found the same points as Chord-with-PGL-Dijkstra:\n";
      // read the points Chord found (from the output file)
      filename.str(std::string());
      filename << "results/query-" << sourceNodeId << "-" << targetNodeId << "-Chord.txt";
      std::vector< pa::PointAndSolution<evns::Path> > rcd = evns::readPointsFromFile(filename.str());
      if ( (rcach.size() == rcd.size()) && std::equal(rcach.begin(), rcach.end(), rcd.begin()) ) 
        std::cout << "  true\n";
      else 
        std::cout << "  false\n";
      // write points to file
      filename.str(std::string());
      filename << "results/query-" << sourceNodeId << "-" << targetNodeId << "-Chord-AStar.txt";
      std::cout << "  Printing points to file " << filename.str() << " ..." << std::endl;
      evns::printPointsToFile(rcach.begin(), rcach.end(), filename.str());
    }

    // The child will now exit. Another child will be made for the next query.
    exit(0);
  }
  else {
    // parent process
    int childStatus;
    pid_t rpid;
    setitimer(ITIMER_REAL, &timeout, NULL);
    rpid = waitpid(pid, &childStatus, 0);
    if (rpid > 0) {
      // child process successfully terminated on its own:
      // first disable the internal timer
      struct itimerval zeroTimeout;
      zeroTimeout.it_value.tv_sec = 0;
      zeroTimeout.it_value.tv_usec = 0;
      zeroTimeout.it_interval.tv_sec = 0;
      zeroTimeout.it_interval.tv_usec = 0;
      setitimer(ITIMER_REAL, &zeroTimeout, NULL);

      // then collect the child's exit status and continue to the next iteration
      waitpid(pid, &childStatus, 0);
      return;
    }
    else if ( (rpid == -1) && (errno == EINTR) ) {
      // waitpid stopped because of an interrupt (we assume it was the alarm):
      // first kill the child process 
      kill(pid, SIGTERM);
      waitpid(pid, &childStatus, 0);
      
      // then print a message and continue to the next iteration
      std::cout << "Killed it. Over 1 hour.\n" << std::endl;
    }
    else {
      // we should never get here
      std::cerr << "Unexpected error during waitpid()... Exiting" << std::endl;
      exit(1);
    }
  }
}


//! The experiments' main function.
//!
int 
main(int argc, char * argv[])
{
  // ----- Defaults -----
  std::string mapPath = "/home/nitsas/Programming/workspace/diplomatiki/dimacs-challenge-graphs/DIMACS9/";
  std::string mapName = "NY";
  std::string coordinatesFilename, graphFilename, distanceFilename, travelTimeFilename;
  std::string queriesFilePath = "/home/nitsas/Programming/workspace/diplomatiki/pareto-approximator/experiments/vs_namoa_star/";
  std::string queriesFilename = "queries.txt";
  bool usingDimacs10Graph = false;

  // ----- Get the map name plus create the full path strings and print them -----

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
    mapPath = "/home/nitsas/Programming/workspace/diplomatiki/dimacs-challenge-graphs/DIMACS10/";
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

  // ----- Parse the queries file -----

  // print the name of the file containing the queries
  std::cout << "Queries file: \n" << (queriesFilePath + queriesFilename) << std::endl;

  // read (parse) the queries file
  std::list< std::pair<unsigned int, unsigned int> > queries;
  queries = evns::parseQueriesFile(queriesFilePath + queriesFilename);

  // print the queries we just read
  std::cout << "Read " << queries.size() 
            << (queries.size() == 1 ? " query:\n" : " queries: \n");
  std::list< std::pair<unsigned int, unsigned int> >::const_iterator qi;
  for (qi = queries.begin(); qi != queries.end(); ++qi) {
    std::cout << qi->first << " " << qi->second << "\n";
  }
  std::cout << std::endl;

  // ----- Make the problem instance (and read the graph) -----

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

  // ----- Make a results directory (if one does not exist) -----
  // Make a directory where all the files containing points of Pareto 
  // and convex Pareto sets will be stored.
  mkdir("results", S_IRWXU | S_IRGRP | S_IROTH | S_IXGRP | S_IXOTH);

  // ----- Stuff needed for the alarm and signal handling. -----

  struct sigaction sa;
  sa.sa_handler = sigalrmHandler;  // handler
  sa.sa_flags = 0;
  sigemptyset(&sa.sa_mask);

  // set the handler for SIGALRM
  if (sigaction(SIGALRM, &sa, NULL) == -1) {
    char * errorMessage = strerror(errno);
    std::cerr << (errorMessage ? errorMessage : "") << "\n";
    exit(1);
  }

  // make an itimerval struct - will use it with setitimer()
  struct itimerval timeout;
  // set the timeout to 1 hour, do not reset automatically after expiration
  // - timeout.it_interval is the next interval to set the timer to, 
  //   after it expires the first time - we set it to 0 so that the timer 
  //   won't reset automatically
  timeout.it_value.tv_sec = 3600;
  timeout.it_value.tv_usec = 0;
  timeout.it_interval.tv_sec = 0;
  timeout.it_interval.tv_usec = 0;

  // ----- Run queries using Chord ----- 

  unsigned int numObjectives = 3;        // CHANGE-HERE
  bool useNamoaStar = false, useAStar = false;

  std::cout << "\nRunning queries (" << mapName << " map, " << numObjectives << " objectives):\n";
  std::cout << "----------------------------------------------------" << std::endl;

  for (qi = queries.begin(); qi != queries.end(); ++qi) {
    std::cout << "FROM NODE " << qi->first << " TO NODE " << qi->second << std::endl;

    // CHANGE--HERE
    // Run a query using the Chord algorithm with PGL's Dijkstra implementation.
    useNamoaStar = false;
    useAStar = false;
    forkAndRunQuery(problem, qi->first, qi->second, timeout, numObjectives, useNamoaStar, useAStar);

    // clean node attributes
    problem.cleanNodeAttributes();

    // CHANGE--HERE
    // Run a query using the Chord algorithm with PGL's A* implementation.
    useNamoaStar = false;
    useAStar = true;
    forkAndRunQuery(problem, qi->first, qi->second, timeout, numObjectives, useNamoaStar, useAStar);

    // clean node attributes
    problem.cleanNodeAttributes();

    std::cout << std::endl;
  }

  // ----- Run the same queries using NAMOA* ----- 

  std::cout << "\nWill now run the same queries using PGL's NAMOA*.\n";
  std::cout << "----------------------------------------------------" << std::endl;

  for (qi = queries.begin(); qi != queries.end(); ++qi) {
    std::cout << "FROM NODE " << qi->first << " TO NODE " << qi->second << std::endl;

    useNamoaStar = true;
    forkAndRunQuery(problem, qi->first, qi->second, timeout, numObjectives, useNamoaStar);

    // clean node attributes
    problem.cleanNodeAttributes();

    std::cout << std::endl;
  }

  return 0;
}


/*! 
 *  @}
 */
