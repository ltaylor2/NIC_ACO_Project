#include "ACO.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>

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
    
    if (argc != 5) {
        std::cout << "USAGE fileName algorithm [ACS | EAS] alpha beta" << std::endl;
        return -1;
    }
    
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
    
    std::string fileName(argv[1]);
    double alpha = atof(argv[3]);
    double beta = atof(argv[4]);

    std::ofstream file;
    file.open("ACO_tests.csv", std::ofstream::out);

    for (int i = 0; i < 3; i++) {
        ACO aco(numIterations, numAnts, alpha, beta, rho, epsilon, q,
                elitismFactor, alg, endCondition, endPercent, maxTime, fileName);
        file << algorithm << "," << alpha << "," << beta << "," << aco.getBestTourPercentage() << std::endl;
        std::cout << "ACO run finished! Best Tour Percentage: " << aco.getBestTourPercentage() << std::endl;
    }

    file.close();

    
    return 0;
}