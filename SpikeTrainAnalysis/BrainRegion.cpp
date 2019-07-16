#include "BrainRegion.h"
#include <string>

BrainRegion::BrainRegion(std::string FileName,int Offset)
	:
	DataFile(FileName, std::ios::binary)
{
	Units.reserve(GetUnitNumber(Offset));
}

int BrainRegion::GetUnitNumber(int Offset)
{
	unsigned short UnitNumber;
	DataFile.seekg(Offset, DataFile.beg);
	char * dataPtr = reinterpret_cast<char*>(&UnitNumber);
	DataFile.read(dataPtr, 2);

	return UnitNumber;
}
