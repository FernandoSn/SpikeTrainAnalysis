#include "BrainRegion.h"
#include <string>

BrainRegion::BrainRegion(std::ifstream* DataFile, unsigned short UnitNumber, int SizePos, int TrainPos)
	:
	Units(UnitNumber)
{
	unsigned int SizeData; // var for getting each param out of the DataFile
	char * SizeDataPtr = reinterpret_cast<char*>(&SizeData);
	DataFile->seekg(SizePos, DataFile->beg);

	for (auto i = Units.begin(), end = Units.end(); i<end ; ++i)
	{
		DataFile->read(SizeDataPtr, 4);
		i->resize(SizeData);
	}

	//double TrainData; // var for getting each param out of the DataFile
	//char * TrainDataPtr = reinterpret_cast<char*>(&TrainData);
	DataFile->seekg(TrainPos, DataFile->beg);

	for (auto i = Units.begin(), end = Units.end(); i<end; ++i)
	{
		DataFile->read(reinterpret_cast<char*>(i->data()), (long long)i->size() * 4);
	}
}

std::vector<std::vector<uint32_t>>& BrainRegion::RUnits()
{
	return Units;
}

std::vector<std::vector<uint32_t>> BrainRegion::UnitsCopy()
{
	return Units;
}
