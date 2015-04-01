#include "ACO.h"

int main(int argc, char** argv)
{
	ACO aco(1000, 1, 3, 0.1, 0.1, 0.04, 0.9, AlgType::elitism, EndCondition::maxIterations, 1, "u2152.tsp");
	return 0;
}