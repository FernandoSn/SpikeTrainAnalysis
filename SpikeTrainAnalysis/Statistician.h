#pragma once

#include "Experiment.h"
#include "BrainRegion.h"

class Statistician
{
public:

	Statistician(std::string FileName);
	Statistician(std::string FileName, double Interval);

private:

	void SetStimLockedSpikes();
	void SetPREXLockedSpikes(double Interval);
	void SpikeTrainCorrelation(const std::vector<double>& reference, const std::vector<double>& target, int BinSize, int epoch);

private:

	Experiment OdorEx;
	BrainRegion Reference;
	BrainRegion Target;

	std::vector<std::vector<double>> StimLockedSpikesRef;
	std::vector<std::vector<double>> StimLockedSpikesTar;
};