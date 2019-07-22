#pragma once

#include "Experiment.h"
#include "BrainRegion.h"
#include <random>
#include <atomic>
#include <mutex>

class Statistician
{
public:

	Statistician(std::string FileName, int BinSize, int Epoch);
	Statistician(std::string FileName, double Interval, int BinSize, int Epoch);
	void RunThreadPool(int ResampledSets, unsigned char ResamplingMethod, double ZThresh, bool ExcZeroLag);
	void RunSingleThread(int ResampledSets, unsigned char ResamplingMethod, double ZThresh, bool ExcZeroLag);

private:

	void SetStimLockedSpikes();
	void SetPREXLockedSpikes(double Interval);
	void SpikeTrainCorr(const std::vector<double>& reference, const std::vector<double>& target, std::vector<unsigned int>& Spikes, unsigned int& Count);
	void SpikeTrainJitter(const std::vector<double>& reference, std::vector<double> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count);
	void SpikeTrainShuffle(const std::vector<double>& reference, std::vector<double> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count);
	void MasterSpikeCrossCorrWorker(int Stimulus, int ResampledSets, unsigned char ResamplingMethod, double ZThresh, bool ExcZeroLag);
	void MasterSpikeCrossCorrDeprecated(int ResampledSets, unsigned char ResamplingMethod, double ZThresh, bool ExcZeroLag);

	//void SpikeTrainShift(); // I dont know if Im gonna implement shift, seems that is not very useful for my actual experiment.

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

	std::random_device Rd;
	std::default_random_engine Generator;

	std::mutex mu;

};

constexpr unsigned char SHUFFLING = 0;
constexpr unsigned char JITTERING = 1;