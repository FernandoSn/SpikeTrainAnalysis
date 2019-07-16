#pragma once

#include "Experiment.h"
#include "BrainRegion.h"

class Statistician
{
public:

	Statistician(std::string FileName);

private:

	void SetStimLockedSpikes(char RegionType);

private:

	Experiment OdorEx;
	BrainRegion Reference;
	BrainRegion Target;

	std::vector<std::vector<double>> StimLockedSpikesRef;
	std::vector<std::vector<double>> StimLockedSpikesTar;
};