#include <string>
#include <utility>

class Graph {
public:
    Graph(std::string filename);
    ~Graph();

    int getNumNodes() const { return numNodes; }
    double getWeight(int i, int j) const { return table[i][j].first; }
    double getPheromone(int i, int j) const { return table[i][j].second; }
    void setPheromone(int i, int j, double pheromone) { table[i][j].second = pheromone; }

private:
    int numNodes;
    std::pair<double, double>** table;
};
