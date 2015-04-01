#pragma once

class Cities {

extern std::string filepath;

public:
	Cities(std::string filepath_);
	int getNumCities();
	double getDistance(int city1, int city2);
	void printCities();

private:
	void readCities();
	std::vector<Locations> cityList;
};

struct Locations {
	double xCoord;
	double yCoord;
};
