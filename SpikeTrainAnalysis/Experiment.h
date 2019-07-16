#pragma once
#include <fstream>

class Experiment
{
public:

	Experiment(std::string FileName);
	unsigned short GetStimuli();
	unsigned short GetMagnitudes();
	unsigned short GetTrials();
	unsigned short GetUnitsRef();
	unsigned short GetUnitsTar();

	int GetRefSizePos();
	int GetRefTrainPos();
	int GetTarSizePos();
	int GetTarTrainPos();

	int GetTimesOnPos();
	int GetTimesOffPos();
	int GetPREXTimesPos();

	std::ifstream* RDataFile();

private:

	std::ifstream DataFile;

	unsigned short Stimuli;
	unsigned short Magnitudes;
	unsigned short Trials;
	unsigned short UnitsRef;
	unsigned short UnitsTar;

	int RefSizePos;
	int RefTrainPos;
	int TarSizePos;
	int TarTrainPos;

	int TimesOnPos;
	int TimesOffPos;
	int PREXTimesPos;

};
