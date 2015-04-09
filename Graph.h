#pragma once

#include <string>
#include <utility>
#include <unordered_map>

/*
 
 ****************************
 *          Graph           *
 ****************************
 Daniel Cohen, Josh Imhoff, and Liam Taylor. 2015. CS3445, Bowdoin College.
 
*/

// Graph class for representing a TSP problem
// Stores weight (distance from city to city) and pheromone in a two dimensional
// matrix indexed by city number, allows O(1) access to both weight and pheromone
class Graph {
public:
    // Constructor
    // @param filename -- the name of the TSP to load
    Graph(std::string filename);
    ~Graph();
    
    // Getters
    // O(1) access
    int getNumNodes() const { return numNodes; }
    inline double getWeight(int i, int j) const { return table[i][j].first; }
    inline double getPheromone(int i, int j) const { return table[i][j].second; }
    double getOptimum() const { return optimalTourWeight; }

    // Setters
    // O(1) access
    void setPheromone(int i, int j, double pheromone) { table[i][j].second = pheromone; }

private:
    int numNodes;
    double optimalTourWeight;
    std::pair<double, double>** table; // NOTE pair.first = weight, pair.second = pheromone

    // NOTE hardcoded based on printout
    static std::unordered_map<std::string, double> optimalTours;
};
