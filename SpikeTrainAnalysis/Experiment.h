#pragma once
#include <fstream>
#include <vector>

class Experiment
{
public:

	Experiment(std::string FileName, bool IsSpontaneous);
	~Experiment();

	int GetStimuli();
	int GetMagnitudes();
	int GetTrials();
	int GetUnitsRef();
	int GetUnitsTar();

	int GetRefSizePos();
	int GetRefTrainPos();
	int GetTarSizePos();
	int GetTarTrainPos();

	int GetTimesOnPos();
	int GetTimesOffPos();
	int GetPREXTimesPos();

	std::ifstream* RDataFile();
	std::vector<double>& GetStimOn();
	std::vector<double>& GetStimOff();
	std::vector<double>& GetPREXTimes();

private:

	void SetNumericalParams(bool IsSpontaneous);
	void SetExpDataVectors(bool IsSpontaneous);

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

	std::vector<double> StimOn;
	std::vector<double> StimOff;
	std::vector<double> PREXTimes;

};
