#include "ACO.h"

#include <vector>
#include <limits>
#include <cmath>
#include <iostream>

ACO::ACO(int numIterations_, double alpha_, double beta_, 
		double rho_, int elitism_, double epsilon_,
		double tao_, double q_, std::string algType_,
		std::string endCondition_, std::string filepath_)
	: numIterations(numIterations_), alpha(alpha_), beta(beta_),
	  rho(rho_), elitism(elitism_), epsilon(epsilon_),
	  tao(tao_), q(q_), algType(algType_), endCondition(endCondition_),
	  graph(filepath_)
{
	run();
}

void ACO::run() {
	int numAnts = graph.getNumNodes();

	graph.initPheromone();

	std::vector<std::vector<int>> ants(numAnts, std::vector<int>(graph.getNumNodes(), 0));
	// initialize all the pheromone levels to 0

	// run through the algorithm
	for (int n = 0; n < numIterations; n++) {
		// evey ant for every iteration
		double bestTourWeight = std::numeric_limits<double>::max();
		int bestAnt = 0;

		for (int i = 0; i < numAnts; i++) {
			double currTourWeight = 0;
			// fill the new unvisited node
			notVisited.clear();
			for (int k = 0; k < graph.getNumNodes(); k++) {
				notVisited.insert(k);
			}

			// pick a first city and mark it as visited
			ants[i][0] = rand() % graph.getNumNodes();
			notVisited.erase(ants[i][0]);

			// assemble the ant's tour under ACO rules
			for (int j = 0; j < graph.getNumNodes() - 1; j++) {

				bool destPicked = false;

				// find the min/max city, and pick it with prob q;
				double random = static_cast<double>(rand()) / RAND_MAX;
				if (random < q) {
					int minMaxNode = getMinMaxNode(ants[i][j]);
					ants[i][j + 1] = minMaxNode;
					notVisited.erase(minMaxNode);
					destPicked = true;
				}
				if (!destPicked) {
					// otherwise, find the next city according to ant system
					double denominator = getProbDenominator(ants[i][j]);
					double missed = 0.0;
					random = static_cast<double>(rand()) / RAND_MAX;

					// pick the next city w/ prob + pheromone
					for (auto it = notVisited.begin(); it != notVisited.end(); it++) {
						double prob = getProbNumerator(ants[i][j], *it)/ denominator;
						
						if (random <= prob + missed) {
							ants[i][j + 1] = *it;
							notVisited.erase(*it);
							destPicked = true;
							break;
						}
						missed += prob;
					}	
				}
				if (!destPicked) {
					ants[i][j + 1] = *(notVisited.rbegin());
				}

				currTourWeight += graph.getWeight(ants[i][j], ants[i][j + 1]);

				// now that we've picked the edge for this section of the tour,
				// wear away pheromone according to ACO rules
				double acoPheromoneLevel = (graph.getPheromone(ants[i][j], ants[i][j+1]) * (1-epsilon)) + (epsilon * tao);
				if (acoPheromoneLevel < 0.0)
					acoPheromoneLevel = 0.0;
				else if (acoPheromoneLevel > 1.0)
					acoPheromoneLevel = 1.0;

				graph.setPheromone(ants[i][j], ants[i][j + 1], acoPheromoneLevel);

			//	std::cout << "From: " << ants[i][j] << " to: " << ants[i][j+1] << std::endl;
			}

			// set the last leg of the tour
			double acoPheromoneLevel = (graph.getPheromone(ants[i][graph.getNumNodes() - 1], ants[i][0]) * (1-epsilon)) + (epsilon * tao);
			if (acoPheromoneLevel < 0.0)
				acoPheromoneLevel = 0.0;
			else if (acoPheromoneLevel > 1.0)
					acoPheromoneLevel = 1.0;
			
			graph.setPheromone(ants[i][graph.getNumNodes() - 1], ants[i][0], acoPheromoneLevel);

			currTourWeight += graph.getWeight(ants[i][graph.getNumNodes() - 1], ants[i][0]);
			
			//std::cout << "CurrTour:" << currTourWeight << std::endl;

			if (currTourWeight < bestTourWeight) {
				bestTourWeight = currTourWeight;
				bestAnt = i;
			}
		}

		// now that all the tours are assembled, plant new pheromone levels according to ACO rules
		// for all tours, change by (1 - rho)*tao
		for (int c = 0; c < graph.getNumNodes(); c++) {
			for (int d = 0; d < graph.getNumNodes(); d++) {
				graph.setPheromone(c, d, (1-rho) * graph.getPheromone(c, d));
			}
		}

		// then add pheromone to edges found in the best tour
		for (int l = 0; l < graph.getNumNodes(); l++) {
			int dest = l;
			if (l == graph.getNumNodes())
				dest = ants[bestAnt][0];
			graph.setPheromone(ants[bestAnt][l], ants[bestAnt][dest], 
							   graph.getPheromone(ants[bestAnt][l], ants[bestAnt][dest]) + (rho * (1 / bestTourWeight)));
		}

		std::cout << "Iteration: " << n << "  Best Tour Length: " << bestTourWeight << std::endl;

	}
}

double ACO::getProbDenominator(int curr) {
	double sum = 0.0;
	//std::cout << "# cities remaining: " << notVisited.size() << std::endl;
	for (auto it = notVisited.begin(); it != notVisited.end(); it++) {
		sum += pow(graph.getPheromone(curr, *it), alpha)*pow(graph.getWeight(curr, *it), beta);
	}

	return sum;
}

double ACO::getProbNumerator(int curr, int dest)
{
	return pow(graph.getPheromone(curr, dest), alpha) * 
	       pow(graph.getWeight(curr, dest), beta);
}

int ACO::getMinMaxNode(int curr)
{
	int minMaxDest = 0;
	int minMaxFactor = 0;
	for (auto it = notVisited.begin(); it != notVisited.end(); it++) {
		int factor = graph.getPheromone(curr, *it) *
					 pow(graph.getWeight(curr, *it), beta);

		if (factor > minMaxFactor) {
			minMaxDest = *it;
			minMaxFactor = factor;
		}
	}

	return minMaxDest;
}
