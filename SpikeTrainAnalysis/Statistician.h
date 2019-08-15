#pragma once

#include "Experiment.h"
#include "BrainRegion.h"
#include <random>
#include <atomic>
#include <mutex>

class Statistician
{
public:

	Statistician(std::string FileName, int BinSize, int Epoch, bool IsSpontaneous);
	Statistician(std::string FileName, int BinSize, int Epoch, double Interval);
	void RunThreadPool(int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double ZThresh, bool ExcZeroLag);
	void RunSingleThread(int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double ZThresh, bool ExcZeroLag);

private:

	void SetStimLockedSpikes();
	void SetPREXLockedSpikes(double Interval);
	void SpikeTrainCorr(const std::vector<double>& reference, const std::vector<double>& target, std::vector<unsigned int>& Spikes, unsigned int& Count);
	void SpikeTrainJitter(const std::vector<double>& reference, const std::vector<double>& target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count);
	void SpikeTrainShuffle(const std::vector<double>& reference, std::vector<double> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count);
	void MasterSpikeCrossCorrWorker(int Stimulus, int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double ZorPVal, bool ExcZeroLag);
	
	template <typename T>
	void WriteToFileWorkerT(std::ofstream& CorrFile, std::vector<T>& CorrVec)
	{
		std::for_each(CorrVec.begin(), CorrVec.end(),
			[&CorrFile](T& Bin)
			{
				CorrFile << Bin << ", ";
			});
		CorrFile << 1 << ", ";
	}

	template <class Iterator, class DataType>
	void STALowerBoundT(Iterator& FI, const Iterator& LI, DataType Limit)
	{
		//FI is the start iterator.
		//LI should always be the end It of the container.

		while (FI < LI && *FI < Limit)
		{
			FI += 1;
		}
	}

	template <class Iterator, class DataType>
	void STAUpperBoundT(Iterator& FI, const Iterator& LI, DataType Limit)
	{
		while (FI < LI && *FI <= Limit)
		{
			FI += 1;
		}
	}

private:

	Experiment OdorEx;
	BrainRegion Reference;
	BrainRegion Target;

	std::vector<std::vector<double>> StimLockedSpikesRef;
	std::vector<std::vector<double>> StimLockedSpikesTar;

	std::atomic<int> BinSize;
	std::atomic<int> Epoch;
	std::atomic<int> NoBins;
	std::atomic<double> BinSizeSec;
	std::atomic<double> EpochSec;

	std::random_device Rd;
	std::default_random_engine Generator;

	std::mutex mu;

};

constexpr unsigned char SHUFFLING = 0;
constexpr unsigned char JITTERING = 1;

constexpr unsigned char ZTEST = 0;
constexpr unsigned char PERMUTATIONTEST = 1;

#define STALowerBoundTExc(FI, LI, Limit) STAUpperBoundT(FI, LI, Limit);
#define STAUpperBoundTExc(FI, LI, Limit) STALowerBoundT(FI, LI, Limit);