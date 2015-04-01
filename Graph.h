#include <string>
#include <utility>
#include <unordered_map>

class Graph {
public:
    Graph(std::string filename);
    ~Graph();

    int getNumNodes() const { return numNodes; }
    double getWeight(int i, int j) const { return table[i][j].first; }
    double getPheromone(int i, int j) const { return table[i][j].second; }
    void setPheromone(int i, int j, double pheromone) { table[i][j].second = pheromone; }
    double getOptimum() const { return optimalTourWeight; }

private:
    int numNodes;
    double optimalTourWeight;
    std::pair<double, double>** table;

    static std::unordered_map<std::string, double> optimalTours;
};
