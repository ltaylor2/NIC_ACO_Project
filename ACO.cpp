#include "ACO.h"

#include <vector>
#include <limits>
#include <cmath>
#include <iostream>
#include <ctime>
#include <time.h>

/*
 
 ****************************
 *            ACO           *
 ****************************
 Daniel Cohen, Josh Imhoff, and Liam Taylor. 2015. CS3445, Bowdoin College.
 
*/

ACO::ACO(int numIterations_, int numAnts_,
        double alpha_, double beta_,
        double rho_, double epsilon_,
        double q_, int elitismFactor_,
        AlgType algType_, EndCondition endCondition_,
        double endPercent_, double maxTime_,
        std::string filepath_)
    : numIterations(numIterations_), numAnts(numAnts_), 
      alpha(alpha_), beta(beta_),
      rho(rho_), epsilon(epsilon_),
      q(q_), elitismFactor(elitismFactor_),
      algType(algType_), endCondition(endCondition_),
      endPercent(endPercent_), maxTime(maxTime_),
      graph(filepath_)
{
    run();
}

void ACO::run() {
    // better start timing!
    srand(time(NULL));
    std::clock_t start = std::clock();

    tao = 0;

    if (algType == AlgType::acs) {
        // get nearest neighbor tour and calculate tao
        double nnTourWeight = getTourLengthNN();
        tao = 1/(numAnts * nnTourWeight);
        
        // initialize pheremone levels
        for (int i = 0; i < graph.getNumNodes(); i++) {
            for (int j = 0; j < graph.getNumNodes(); j++) {
                if (i != j) {
                    graph.setPheromone (i, j, tao);
                }
            }
        }

    }    
    
    // ant is a std::vector<int> that stores cities, representing a tour
    std::vector<std::vector<int>> ants(numAnts, std::vector<int>(graph.getNumNodes(), 0));
    
    // best tour is a std::vector<int> that stores cities, best tour found so far
    std::vector<int> bestTour(graph.getNumNodes(), 0);
    
    // length of bestTour
    bestTourWeight = std::numeric_limits<double>::max();

    // optimal tour as recorded in document and saved in Graph.cpp
    double optimum = graph.getOptimum();
    
    // tourWeights is a vector of weights of current tour per ant used in Elitist ACO
    std::vector<double> tourWeights(numAnts, 0);

    int n = 0;
    while (n < numIterations) {
        
        // NOTE only update n if stopping condition is maxIterations or all
        //      otherwise no reason to keep track of iterations
        if (endCondition == EndCondition::maxIterations || endCondition == EndCondition::all)
            n++;
        
        // loop over ants
        for (int i = 0; i < numAnts; i++) {
            double currTourWeight = 0;
            
			// no city has been visited as of now
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
                        // being picked NOTE linear time
                        if (random <= prob + missed) {
                            ants[i][j + 1] = *it;
                            destPicked = true;
                            break;
                        }
                        missed += prob;
                    }
                }
                
                // if none have been picked, pick the last city in the list
                // NOTE would be result of janky floating point math
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
                    bestTour[b] = ants[i][b];
            }
            
            if (endCondition == EndCondition::afterPercentageOfOptimal || endCondition == EndCondition::all) {
                if (optimum != 0 && bestTourWeight * endPercent <= optimum)
                    break;
            }
        }

        // now that all the tours are assembled
        allBestTourWeights.push_back(bestTourWeight / optimum);

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

        std::clock_t end = std::clock();
        double currTime = (end - start) / (double) CLOCKS_PER_SEC;
        
        // print status update
        std::cout << "Iteration: " << n << "  Best Tour Percentage: " << bestTourWeight / optimum << " Time: " << currTime << std::endl;
        
        // check the current conditions, and end the run if they match
        // the current end condition specifications
        if (endCondition == EndCondition::afterPercentageOfOptimal || endCondition == EndCondition::all) {
            if (bestTourWeight / optimum <= endPercent) {
                std::cout << "Percent of optimum achieved." << std::endl;
                break;
            }
        }

        if (endCondition == EndCondition::time || endCondition == EndCondition::all) {
            if (currTime >= maxTime) {
                std::cout << "Max time completed." << std::endl;
                break;
            }
        }
    }
}

// denom only needs be calculated once for ant
// @param curr -- refers to ant
// NOTE linear time in size of notVisited set
double ACO::getProbDenominator(int curr) {
    double sum = 0.0;
    //std::cout << "# cities remaining: " << notVisited.size() << std::endl;
    for (auto it = notVisited.begin(); it != notVisited.end(); it++) {
        sum += pow(graph.getPheromone(curr, *it), alpha)*pow(1 / graph.getWeight(curr, *it), beta);
    }

    return sum;
}

// num needs to be calculated once per city in tour
// @param curr -- refers to ant
// @param dest -- refers to city being considered
// NOTE O(1)
double ACO::getProbNumerator(int curr, int dest)
{
    return pow(graph.getPheromone(curr, dest), alpha) * 
           pow(1 / graph.getWeight(curr, dest), beta);
}

// get minMaxNode for ACS
// @param curr -- refers to ant
// NOTE linear time in size of notVisited set
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

// get the nearest neighbor tour to find value
// of tao in ACS variant
double ACO::getTourLengthNN()
{
    double tourLength = 0;
    int curNode = (int)rand() % graph.getNumNodes();
    std::vector<int> nearestNeighborTour;
    //nearestNeighborTour.push_back(curNode);
    
    std::set<int> notVisitedNN;
    for (int k = 0; k < graph.getNumNodes(); k++) {
        notVisitedNN.insert(k);
    }
    notVisitedNN.erase(curNode);
    
    while (!notVisitedNN.empty()) {
        int closestNode = std::numeric_limits<int>::max();
        double closestWeight = std::numeric_limits<double>::max();
        for (auto it = notVisitedNN.begin(); it != notVisitedNN.end(); it++) {
            double curWeight = graph.getWeight(curNode, *it);
            if (curWeight < closestWeight) {
                closestWeight = curWeight;
                closestNode = *it;
            }
        }
        tourLength += closestWeight;
        notVisitedNN.erase(closestNode);
        curNode = closestNode;
    }
    return tourLength;
}