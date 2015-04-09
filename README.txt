**************************************************************************************
Ant Colony Optimization for Traveling Salesman Problem:A Comparison of Ant Colony System and Elitist Ant System and an Exploration of the Alpha and Beta***************************************************************************************
Dan Cohen, Josh Imhoff, and Liam Taylor
April 8, 2015
Project 3
CS3445

This project is designed to test differences between Ant Colony System (ACS) and Elitist Ant System (EAS) variants of Ant Colony Optimization (ACO). Due to our experimental design, all parameters are hard-coded in the main method (Main.cpp), but can be changed there. The only cmd line arguments that change parameters are alpha and beta.

TO RUN:
After running the makefile in the project folder directory, run the program from the command line using the following syntax:

“./aco”

then,
	* fileName: .tsp file directory/name

	* algorithm (decides what ACO variant to use):
				“ACS” -  Ant Colony System
				“EAS” -  Elitist Ant System
		
	* alpha: effectively the weight of pheromone influence on tour assembly

	
	* beta: effectively the weight of edge-heuristic information (1/length) on tour assembly

For example, to run the ACS variant on an example .tsp file with alpha = 1.5 and beta = 5, enter:

“./aco u2152.tsp ACS 1.5 5”
