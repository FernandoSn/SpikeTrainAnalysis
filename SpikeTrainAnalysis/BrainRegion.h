#pragma once

#include <vector>
#include <fstream>

class BrainRegion
{
public:

	BrainRegion(std::ifstream* DataFile, unsigned short UnitNumber, int SizePos, int TrainPos);
	std::vector<std::vector<double>>& RUnits();

private:

	
private:

	std::vector<std::vector<double>> Units;
};
