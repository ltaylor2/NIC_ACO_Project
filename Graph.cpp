#include "Graph.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

Graph::Graph(std::string filename)
{ 
    std::fstream file(filename);
    std::stringstream ss;

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
            } else if (line.compare("NODE_COORD_SECTION") == 0 ||
                       line.compare("DISPLAY_DATA_SECTION") == 0) {
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
void Graph::initPheromone()
{
    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) {
            setPheromone(i, j, 1);
        }
    }
}

Graph::~Graph()
{
    for (int i = 0; i < numNodes; i++)
        delete[] table[i];
    delete[] table;
}
