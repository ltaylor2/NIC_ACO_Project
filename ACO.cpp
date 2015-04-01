#include "ACO.h"

#include <vector>
#include <limits>
#include <cmath>
#include <iostream>

ACO::ACO(int numIterations_, double alpha_, double beta_, 
		double rho_, double epsilon_,
		double tao_, double q_, AlgType algType_,
		EndCondition endCondition_, double endPercent_, std::string filepath_)
	: numIterations(numIterations_), alpha(alpha_), beta(beta_),
	  rho(rho_), epsilon(epsilon_),
	  tao(tao_), q(q_), algType(algType_), endCondition(endCondition_),
      endPercent(endPercent_), graph(filepath_)
{
	run();
}

void ACO::run() {
    
    // number of ants is always the number of nodes
	//int numAnts = graph.getNumNodes();
    int numAnts = 20;
    
	// elitism factor is always set as the number of ants
	int elitismFactor = numAnts;

	std::vector<std::vector<int>> ants(numAnts, std::vector<int>(graph.getNumNodes(), 0));
	std::vector<int> bestTour(graph.getNumNodes(), 0);
	double bestTourWeight = std::numeric_limits<double>::max();

    double optimum = graph.getOptimum();
    
	// run through the algorithm
    //keep track of all the ants tour weights for elistist updates from each ant
    std::vector<double> tourWeights(numAnts, 0);

    int n = 0;
	while (n < numIterations) {
        if (endCondition == EndCondition::maxIterations || endCondition == EndCondition::all)
            n++;
            
		for (int i = 0; i < numAnts; i++) {
			double currTourWeight = 0;
            
			// fill the new unvisited node
			notVisited.clear();
			for (int k = 0; k < graph.getNumNodes(); k++) {
				notVisited.insert(k);
			}

			// pick a first city and mark it as visited
			ants[i][0] = (int)rand() % graph.getNumNodes();

			notVisited.erase(ants[i][0]);

			// assemble the ant's tour
			for (int j = 0; j < graph.getNumNodes() - 1; j++) {
				bool destPicked = false;

				// find the min/max city, and pick it with prob q;
				double random = static_cast<double>(rand()) / RAND_MAX;
                
                // for ant colony system, choose city that maximizes tao * rho ^ beta with probability q
				if (algType == AlgType::acs) {
                    if (random < q) {
                        int minMaxNode = getMinMaxNode(ants[i][j]);
                        ants[i][j + 1] = minMaxNode;
                        destPicked = true;
                    }
                }

				if (!destPicked) {
					// pick the next city according to the ant system
					double denominator = getProbDenominator(ants[i][j]);
					double missed = 0.0;
					random = static_cast<double>(rand()) / RAND_MAX;

					// pick the next city w/ prob + pheromone
					for (auto it = notVisited.begin(); it != notVisited.end(); it++) {
						double prob = getProbNumerator(ants[i][j], *it)/ denominator;
						
                        // pick cities linearly by adding probability of  previously skipped cities to those
                        // being picked
						if (random <= prob + missed) {
							ants[i][j + 1] = *it;
							destPicked = true;
							break;
						}
						missed += prob;
					}
				}
                
                // if none have been picked from the linear prob. algorithm, pick the last city in the list
                // (should have been prob == 1 anyways
				if (!destPicked) {
					ants[i][j + 1] = *(notVisited.rbegin());
				}
                
                // erase the new city from the not visited list
				notVisited.erase(ants[i][j+1]);

                // keep track of the tour weight
				currTourWeight += graph.getWeight(ants[i][j], ants[i][j + 1]);

                // wear away pheromone according to acs rules
                if (algType == AlgType::acs) {
                    // now that we've picked the edge for this section of the tour,
                    // wear away pheromone according to ACO rules
                    double acoPheromoneLevel = (graph.getPheromone(ants[i][j], ants[i][j+1]) * (1-epsilon)) + (epsilon * tao);
                    graph.setPheromone(ants[i][j], ants[i][j + 1], acoPheromoneLevel);
                }
			}
            
            // wear away pheromone on the last leg according to acs rules
            if (algType == AlgType::acs) {
                double acoPheromoneLevel = (graph.getPheromone(ants[i][graph.getNumNodes() - 1], ants[i][0]) * (1-epsilon)) + (epsilon * tao);
                graph.setPheromone(ants[i][graph.getNumNodes() - 1], ants[i][0], acoPheromoneLevel);
            }
            
			// set the last leg of the tour and add its weights
			currTourWeight += graph.getWeight(ants[i][graph.getNumNodes() - 1], ants[i][0]);
			
            // keep track of all the tour weights for elitism
            if (algType == AlgType::elitism)
                tourWeights[i] = currTourWeight;

            // keep track of the best tour so far
			if (currTourWeight < bestTourWeight) {
				bestTourWeight = currTourWeight;
                for (int b = 0; b < graph.getNumNodes(); b++)
                    bestTour[b] = ants [i][b];
			}
            
            if (endCondition == EndCondition::afterPercentageOfOptimal || endCondition == EndCondition::all) {
                if (optimum != 0 && bestTourWeight * endPercent <= optimum)
                    break;
            }
		}

		// now that all the tours are assembled
        
        // run through all ant tours and update elitism ant summation based on their tour lengths
        if (algType == AlgType::elitism) {
            for (unsigned int i = 0; i < tourWeights.size(); i++) {
                for (int l = 0; l < graph.getNumNodes(); l++) {
                    int dest = l + 1;
                    if (dest == graph.getNumNodes())
                        dest = ants[i][0];
                    graph.setPheromone(ants[i][l], ants[i][dest],
                                       graph.getPheromone(ants[i][l], ants[i][dest]) + (1 / tourWeights[i]));
                }
            }
        }

		// for all tours (in both elitism and acs), change by (1 - rho)*tao
		for (int c = 0; c < graph.getNumNodes(); c++) {
			for (int d = 0; d < graph.getNumNodes(); d++) {
				graph.setPheromone(c, d, (1-rho) * graph.getPheromone(c, d));
			}
		}

		// then add pheromone to edges found in the best tour found so far
		for (int l = 0; l < graph.getNumNodes(); l++) {
			int dest = l + 1;
			if (dest == graph.getNumNodes())
				dest = bestTour[0];
            
            // for acs, use rho to update the best tour pheromone
            if (algType == AlgType::acs) {
                graph.setPheromone(bestTour[l], bestTour[dest],
                                   graph.getPheromone(bestTour[l], bestTour[dest]) + (rho * (1 / bestTourWeight)));
            }
            
            // for elitism, use elitism factor to update best tour pheromone
            if (algType == AlgType::elitism) {
                graph.setPheromone(bestTour[l], bestTour[dest],
                                   graph.getPheromone(bestTour[l], bestTour[dest]) + (elitismFactor * (1 / bestTourWeight)));
		
            }
        }

		std::cout << "Iteration: " << n << "  Best Tour Length: " << bestTourWeight << std::endl;
	}
}

double ACO::getProbDenominator(int curr) {
	double sum = 0.0;
	//std::cout << "# cities remaining: " << notVisited.size() << std::endl;
	for (auto it = notVisited.begin(); it != notVisited.end(); it++) {
		sum += pow(graph.getPheromone(curr, *it), alpha)*pow(1 / graph.getWeight(curr, *it), beta);
	}

	return sum;
}

double ACO::getProbNumerator(int curr, int dest)
{
	return pow(graph.getPheromone(curr, dest), alpha) * 
	       pow(1 / graph.getWeight(curr, dest), beta);
}

int ACO::getMinMaxNode(int curr)
{
	int minMaxDest = 0;
	double minMaxFactor = 0;
	for (auto it = notVisited.begin(); it != notVisited.end(); it++) {
		double factor = graph.getPheromone(curr, *it) *
					 pow(1 / graph.getWeight(curr, *it), beta);

		if (factor > minMaxFactor) {
			//std::cout << "factor = " << factor << std::endl;
			minMaxDest = *it;
			minMaxFactor = factor;
		}
	}

	return minMaxDest;
}
