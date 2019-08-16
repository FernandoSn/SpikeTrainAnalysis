#pragma once

#include <vector>
#include <fstream>

class BrainRegion
{
public:

	BrainRegion(std::ifstream* DataFile, unsigned short UnitNumber, int SizePos, int TrainPos);
	std::vector<std::vector<uint32_t>>& RUnits();
	std::vector<std::vector<uint32_t>> UnitsCopy();

private:

	
private:

	std::vector<std::vector<uint32_t>> Units;
};
