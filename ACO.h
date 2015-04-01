#pragma once

#include "Graph.h"
#include <set>

class ACO {

public:
	ACO(int numIterations_, double alpha_, double beta_, 
		double rho_, double epsilon_,
		double tao_, double q_, std::string algType_,
		std::string endCondition_, std::string filepath_);

private:
	void run();

	int numIterations;
	double alpha;
	double beta;
	double rho;
	int elitism;
	double epsilon;
	double tao;
	double q;
	std::string algType;
	std::string endCondition;
	std::string filepath;
	std::set<int> notVisited;
	Graph graph;

	double getProbDenominator(int curr);
	double getProbNumerator(int curr, int dest);
	int getMinMaxNode(int curr);
};