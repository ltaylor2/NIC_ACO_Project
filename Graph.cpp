#include "Graph.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

std::unordered_map<std::string, double> Graph::optimalTours {
    {"d2103.tsp", 80450},
    {"u2152.tsp", 64253},
    {"u2319.tsp", 234256},
    {"pr2392.tsp", 378032},
    {"pcb3038.tsp", 137694},
    {"fl3795.tsp", 28772},
    {"fnl4461.tsp", 182566},
    {"rl5915.tsp", 565530},
    {"rl5934.tsp", 556045},
    {"pla7397.tsp", 23260728},
    {"rl11849.tsp", 923288},
    {"pla85900.tsp", 142382641}
};

Graph::Graph(std::string filename)
{ 
    std::fstream file(filename);
    std::stringstream ss;

    try {
        optimalTourWeight = Graph::optimalTours.at(filename);         
        std::cout << "Optimal tour weight: " << optimalTourWeight << std::endl;
    } catch (const std::out_of_range& oor) {
        optimalTourWeight = 0;
        std::cout << "Optimal tour weight not recorded." << std::endl;
    }

    if (file.is_open()) {
        std::string line;

        // TODO problems can be represented as weight matrices
        //          - full matrix
        //          - upper triangular portion of matrix
        std::pair<double, double>* coords;
        while (getline(file, line)) {
            ss.clear();
            ss.str("");
            ss << line;
            if (line.substr(0,4).compare("NAME") == 0) {
                int k = line.find(':');
                std::string problemName = line.substr(k+2);
            } else if (line.substr(0,9).compare("DIMENSION") == 0) {
                int k = line.find(':');
                numNodes = std::stod(line.substr(k+2));
                table = new std::pair<double, double>*[numNodes];
                for (int i = 0; i < numNodes; i++)
                    table[i] = new std::pair<double, double>[numNodes];
                coords = new std::pair<double, double>[numNodes];
            } else if (line.substr(0,18).compare("NODE_COORD_SECTION") == 0 ||
                       line.substr(0,20).compare("DISPLAY_DATA_SECTION") == 0) {
                for (int i = 0; i < numNodes; i++) {
                    getline(file, line);
                    ss.clear();
                    ss.str("");
                    ss << line;
                    int index;
                    ss >> index;
                    index--;
                    double xCoord, yCoord;
                    ss >> xCoord;
                    ss >> yCoord;
                    coords[i].first = xCoord;
                    coords[i].second = yCoord;
                    for (int j = 0; j < i; j++) {
                        double xCoordDiff = xCoord - coords[j].first;
                        double yCoordDiff = yCoord - coords[j].second;
                        double weight = sqrt(xCoordDiff*xCoordDiff + yCoordDiff*yCoordDiff);
                        table[j][index].first = weight;
                        table[index][j].first = weight;
                        table[j][index].second = 1;
                        table[index][j].second = 1;
                    }
                }
            }
        }
        delete[] coords;
    }

    // for (int i = 0; i < numNodes; i++) {
    //     for (int j = 0; j < numNodes; j++) {
    //         std::cout << i << " ";
    //         std::cout << j << " ";
    //         std::cout << table[i][j].first << " ";
    //         std::cout << std::endl;
    //     }
    // }
}

Graph::~Graph()
{
    for (int i = 0; i < numNodes; i++)
        delete[] table[i];
    delete[] table;
}
