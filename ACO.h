#pragma once

#include "Graph.h"
#include <set>

enum class AlgType {
    elitism,
    acs
};

enum class EndCondition {
    maxIterations,
    afterPercentageOfOptimal,
    all
};

class ACO {
public:
	ACO(int numIterations_, double alpha_, double beta_, 
		double rho_, double epsilon_,
		double tao_, double q_, AlgType algType_,
		EndCondition endCondition_, double endPercent_, std::string filepath_);

private:
	void run();

	int numIterations;
	double alpha;
	double beta;
	double rho;
	double epsilon;
	double tao;
	double q;
	AlgType algType;
	EndCondition endCondition;
    double endPercent;
	std::string filepath;
	Graph graph;

    std::set<int> notVisited;

	double getProbDenominator(int curr);
	double getProbNumerator(int curr, int dest);
	int getMinMaxNode(int curr);
};