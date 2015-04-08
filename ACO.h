#pragma once

#include "Graph.h"
#include <set>
#include <vector>

enum class AlgType {
    elitism,
    acs
};

enum class EndCondition {
    maxIterations,
    afterPercentageOfOptimal,
    time,
    all
};

class ACO {
public:
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
	void run();

	int numIterations;
    int numAnts;
    
	double alpha;
	double beta;
	double tao;
	double rho;
	double epsilon;
	double q;
    int elitismFactor;
	AlgType algType;
	EndCondition endCondition;
    double endPercent;
    double maxTime;
	std::string filepath;
	Graph graph;

    std::set<int> notVisited;

    double bestTourWeight;

    std::vector<double> allBestTourWeights;

	double getProbDenominator(int curr);
	double getProbNumerator(int curr, int dest);
	int getMinMaxNode(int curr);
	double getTourLengthNN();
};