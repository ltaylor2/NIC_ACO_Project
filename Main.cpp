#include "ACO.h"

#include <stdlib.h>
#include <iostream>

/*
 
 ****************************
 *          Main            *
 ****************************
 This program attempts to fulfill Project 3 for CS3445 at Bowdoin College.
 Daniel Cohen, Josh Imhoff, and Liam Taylor. 2015.
 This file includes the main function.
 
*/

// Commandline options:
// @arg1 fileName (.tsp file) (string)
// @arg2 algorithm [ACS=Ant Colony System | EAS=Elitist Ant System] (string)
// @arg3 alpha (double)
// @arg4 beta (double)
int main(int argc, char** argv)
{
    // the majority of the variable for the ACS algorithm are hard-coded here, as this
    // program was used to analyze trials where the only differences were tsp file problem,
    // algorithm type, and alpha and beta levels.
    int numIterations = 1000;
    int numAnts = 20;
	double rho = 0.1;
	double epsilon = 0.1;
	double q = 0.9;
    int elitismFactor = numAnts;
	EndCondition endCondition = EndCondition::time;
    double endPercent = 1.05;
    double maxTime = 900;
    
    // check cmd line args
    if (argc != 5) {
        std::cout << "USAGE fileName algorithm [ACS | EAS] alpha beta" << std::endl;
        return -1;
    }
    
    // set the correct ACO variant
    AlgType alg;
    std::string algorithm(argv[2]);
    if (algorithm == "ACS")
        alg = AlgType::acs;
    else if (algorithm == "EAS")
        alg = AlgType::elitism;
    else {
        std::cout << "Invalid Ant Colony Algorithm. Closing." << std::endl;
        return -1;
    }
    
    // set the last of the parameters, just to be orderly
    std::string fileName(argv[1]);
    double alpha = atof(argv[3]);
    double beta = atof(argv[4]);

    // call ACO constructor, runs algs
    ACO aco(numIterations, numAnts, alpha, beta, rho, epsilon, q,
            elitismFactor, alg, endCondition, endPercent, maxTime, fileName);

    // report final output
    std::cout << "ACO run finished! Best Tour Percentage: " << aco.getBestTourPercentage() << std::endl;
    
    return 0;
}