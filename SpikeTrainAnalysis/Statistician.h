#pragma once

#include "Experiment.h"
#include "BrainRegion.h"
#include <random>

class Statistician
{
public:

	Statistician(std::string FileName, int BinSize, int Epoch);
	Statistician(std::string FileName, double Interval, int BinSize, int Epoch);

private:

	void SetStimLockedSpikes();
	void SetPREXLockedSpikes(double Interval);
	void SpikeTrainCorr(const std::vector<double>& reference, const std::vector<double>& target, std::vector<unsigned int>& Spikes, unsigned int& Count);
	void SpikeTrainJitter();
	void SpikeTrainShift();
	void SpikeTrainShuffle(const std::vector<double>& reference, std::vector<double> target);
	void MasterSpikeCrossCorr();
	void InitInterns();
	void MasterSpikeCrossCorrWorker(int Stimulus);

private:

	Experiment OdorEx;
	BrainRegion Reference;
	BrainRegion Target;

	std::vector<std::vector<double>> StimLockedSpikesRef;
	std::vector<std::vector<double>> StimLockedSpikesTar;

	int BinSize;
	int Epoch;
	int NoBins;
	double BinSizeSec;
	double EpochSec;
	static constexpr int Shuffles = 100;
	static constexpr int Shifts = 30;
	static constexpr int Jitters = 1000;

	struct
	{
		unsigned int Corr;
		unsigned int Jitter;
		unsigned int Shift;
		unsigned int Shuffle;
	}Counts;

	std::vector<unsigned int> SpikesCountCorr;
	std::vector<unsigned int> SpikesCountShift;
	std::vector<unsigned int> SpikesCountShuffle;
	std::vector<unsigned int> SpikesCountJitter;

	std::random_device Rd;
	std::default_random_engine Generator;

};