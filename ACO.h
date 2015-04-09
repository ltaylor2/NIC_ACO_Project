#pragma once

#include "Graph.h"

#include <set>
#include <vector>

// specifies ACO variant types
// acs = Ant Colony System
// elitism = Elitist Ant System
enum class AlgType {
    elitism,
    acs
};

// specifies end condition for the run method
// maxIterations = run for a certain number of iterations
// afterPercentageOfOptimal = run until a good enough solution is found
// time = run for a certain amount of time
// all = whatever first from all of the above
enum class EndCondition {
    maxIterations,
    afterPercentageOfOptimal,
    time,
    all
};


/*
 
 ****************************
 *            ACO           *
 ****************************
 Daniel Cohen, Josh Imhoff, and Liam Taylor. 2015. CS3445, Bowdoin College.
 
*/

class ACO {
public:
    // Constructor
	ACO(int numIterations_, int numAnts_,
        double alpha_, double beta_,
		double rho_, double epsilon_,
        double q_, int elitismFactor_,
        AlgType algType_, EndCondition endCondition_,
        double endPercent_, double maxTime_,
        std::string filepath_);

	double getBestTourWeight() { return bestTourWeight; }
	double getBestTourPercentage() { return bestTourWeight / graph.getOptimum(); }
	std::vector<double>* getAllBestWeights() { return &allBestTourWeights; }

private:
    // Main ACO run function
	void run();

    // Helper functions
    double getProbDenominator(int curr);
	double getProbNumerator(int curr, int dest);
	int getMinMaxNode(int curr);
	double getTourLengthNN();

    int numIterations;          // max number of iterations if EndCondition == maxIterations

    // Parameters
    int numAnts;
	double alpha;       // proportional effect of pheromone on tour building
	double beta;        // proportional effect of heuristic on tour building
	double tao;         // ACS init pheromone level and deposition, assembled NN tour
	double rho;         // pheromone evap
	double epsilon;     // pheromone evap (ACS)
	double q;           // min-max probability (ACS)
    int elitismFactor;  // elitism deposition factor (EAS)
    
    
	AlgType algType;            // elitism, acs
    
	EndCondition endCondition;  // maxIterations, afterPercentageOfOptimal, all
    double endPercent;          // percent of optimal tour, if reached, halt
    double maxTime;             // time until a timed run will stop
    
    // Graph stores weights (euclidean distances) and pheromone levels
	Graph graph;

    // notVisited set stores cities that have not yet been chosen in construction of tour
    std::set<int> notVisited;

    double bestTourWeight;

    std::vector<double> allBestTourWeights;
};