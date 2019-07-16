#include "Experiment.h"

Experiment::Experiment(std::string FileName)
	:
	DataFile(FileName, std::ios::binary)
{
	SetNumericalParams();
	SetExpDataVectors();
}

unsigned short Experiment::GetStimuli()
{
	return Stimuli;
}

unsigned short Experiment::GetMagnitudes()
{
	return Magnitudes;
}

unsigned short Experiment::GetTrials()
{
	return Trials;
}

unsigned short Experiment::GetUnitsRef()
{
	return UnitsRef;
}

unsigned short Experiment::GetUnitsTar()
{
	return UnitsTar;
}

int Experiment::GetRefSizePos()
{
	return RefSizePos;
}

int Experiment::GetRefTrainPos()
{
	return RefTrainPos;
}

int Experiment::GetTarSizePos()
{
	return TarSizePos;
}

int Experiment::GetTarTrainPos()
{
	return TarTrainPos;
}

int Experiment::GetTimesOnPos()
{
	return TimesOnPos;
}

int Experiment::GetTimesOffPos()
{
	return TimesOffPos;
}

int Experiment::GetPREXTimesPos()
{
	return PREXTimesPos;
}

std::ifstream* Experiment::RDataFile()
{
	return &DataFile;
}

std::vector<double>& Experiment::GetStimOn()
{
	return StimOn;
}

std::vector<double>& Experiment::GetStimOff()
{
	return StimOff;
}

std::vector<double>& Experiment::GetPREXTimes()
{
	return PREXTimes;
}

void Experiment::SetNumericalParams()
{
	//Helper Variables to obtain data.
	unsigned short Data; // var for getting each param out of the DataFile
	char * DataPtr = reinterpret_cast<char*>(&Data);
	unsigned int TempData;
	char * TempPtr = reinterpret_cast<char*>(&TempData);


	DataFile.seekg(0, DataFile.beg); //Go to the beginning of DataFile

	DataFile.read(DataPtr, 2);
	Stimuli = Data; //Get stimuli aka valves or different odors

	DataFile.read(DataPtr, 2);
	Magnitudes = Data; //Get magnitudes aka odor concentrations.

	DataFile.read(DataPtr, 2);
	Trials = Data;

	DataFile.read(DataPtr, 2);
	UnitsRef = Data;

	DataFile.read(DataPtr, 2);
	UnitsTar = Data;

	RefSizePos = (int)DataFile.tellg();

	TarSizePos = RefSizePos + UnitsRef * 4; // 4 for 4 bytes numbers.

	RefTrainPos = TarSizePos + UnitsTar * 4;

	TarTrainPos = RefTrainPos;
	for (int i = 0; i < UnitsRef; i++)
	{
		DataFile.read(TempPtr, 4);
		TarTrainPos += TempData * 8; //8 because we are reading 8 byte doubles now.
	}

	TimesOnPos = TarTrainPos;
	for (int i = 0; i < UnitsTar; i++)
	{
		DataFile.read(TempPtr, 4);
		TimesOnPos += TempData * 8; //8 because we are reading 8 byte doubles now.
	}

	TimesOffPos = TimesOnPos + Stimuli * Magnitudes * Trials * 8;

	PREXTimesPos = TimesOffPos + Stimuli * Magnitudes * Trials * 8;
}

void Experiment::SetExpDataVectors()
{
	StimOn.reserve(Stimuli * Magnitudes * Trials);
	StimOff.reserve(Stimuli * Magnitudes * Trials);
	PREXTimes.reserve(Stimuli * Magnitudes * Trials);


	DataFile.seekg(TimesOnPos, DataFile.beg);
	DataFile.read(reinterpret_cast<char*>(StimOn.data()), StimOn.capacity() * 8);

	DataFile.seekg(TimesOffPos, DataFile.beg);
	DataFile.read(reinterpret_cast<char*>(StimOff.data()), StimOff.capacity() * 8);

	DataFile.seekg(PREXTimesPos, DataFile.beg);
	DataFile.read(reinterpret_cast<char*>(PREXTimes.data()), PREXTimes.capacity() * 8);
}
