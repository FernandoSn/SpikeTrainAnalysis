#pragma once

#include <vector>
#include <fstream>

class BrainRegion
{
public:

	BrainRegion(std::string FileName, int Units);

private:

	int GetUnitNumber(int Offset);
	
private:

	std::vector<std::vector<double>> Units;
	std::ifstream DataFile;
};
