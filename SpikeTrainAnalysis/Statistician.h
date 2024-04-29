#pragma once

#include "Experiment.h"
#include "BrainRegion.h"
#include <random>
#include <atomic>
#include <mutex>
#include <iterator>
#include <numeric>

class Statistician
{
public:

	Statistician(std::string FileName, int BinSize, int Epoch, bool IsSpontaneous);
	Statistician(std::string FileName, int BinSize, int Epoch, uint32_t Interval);
	void RunThreadPool(int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double ZThresh, bool ExcZeroLag);
	void RunSingleThread(int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double ZThresh, bool ExcZeroLag);

private:

	void SetStimLockedSpikes();
	void SetPREXLockedSpikes(uint32_t Interval);
	void SpikeTrainCorr(const std::vector<uint32_t>& reference, const std::vector<uint32_t>& target, std::vector<unsigned int>& Spikes, unsigned int& Count);
	void SpikeTrainIntervalJitter(const std::vector<unsigned int>& Spikes, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count);
	void SpikeTrainBasicJitter(const std::vector<uint32_t>& reference, const std::vector<uint32_t>& target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count);
	void SpikeTrainBasicCuJitter(const std::vector<uint32_t>& reference, std::vector<uint32_t> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count);
	void SpikeTrainShuffle(const std::vector<uint32_t>& reference, std::vector<uint32_t> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count);
	void MasterSpikeCrossCorrWorker(int ThreadNo, int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double ZorPVal, bool ExcZeroLag);
	void MasterSpikeCrossCorrSingleThread(int ThreadNo, int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double ZorPVal, bool ExcZeroLag);
	void CloseFiles();
	
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

	template <class T>
	void GetSignifcantCorr(const std::vector<T>& CountVec, bool* SigArray, const std::pair<T, T>& GlobalBands, const std::vector<T>& LPWBand, const std::vector<T>& UPWBand)
	{

		std::vector<bool> CountSigVec(NoBins);
		auto LeadBeg = CountSigVec.end() - (CountSigVec.size() / 2);
		auto LeadEnd = CountSigVec.end();
		auto LagBeg = CountSigVec.begin();
		auto LagEnd = LeadBeg;
		std::vector<bool>::iterator SigLeadPosition[2];
		std::reverse_iterator<std::vector<bool>::iterator> SigLagPosition[2];
		int maxSigBins = 4; //I used 4 normally.

		//Excitatory Correlations

		std::transform(CountVec.cbegin(), CountVec.cend(), CountSigVec.begin(), //Count the bins that cross the threshold!
			[&GlobalBands](const T& Spike) -> bool { return Spike > GlobalBands.second; });

		
		int LeadCount = std::accumulate(LeadBeg, LeadEnd, 0);
		int LagCount = std::accumulate(LagBeg, LagEnd, 0);

		//if ((LeadCount == 1 || LeadCount == 2) && ((LeadCount - *LeadBeg - *(LeadEnd - 1)) > 0))
		if (LeadCount > 0 && LeadCount < maxSigBins)
		{
			SigArray[0] = true;
			SigLeadPosition[0] = std::find(LeadBeg, LeadEnd, 1);
		}

		//if ((LagCount == 1 || LagCount == 2) && ((LagCount - *LagBeg - *(LagEnd - 1)) > 0))
		if (LagCount > 0 && LagCount < maxSigBins)
		{
			SigArray[1] = true;
			SigLagPosition[0] = std::find(std::make_reverse_iterator(LagEnd), std::make_reverse_iterator(LagBeg + 1), 1);
		}

		//Inhibitory Correlations

		std::transform(CountVec.cbegin(), CountVec.cend(), CountSigVec.begin(), //Count the bins that cross the threshold!
			[&GlobalBands](const T& Spike) -> bool { return Spike < GlobalBands.first; });


		LeadCount = std::accumulate(LeadBeg, LeadEnd, 0);
		LagCount = std::accumulate(LagBeg, LagEnd, 0);

		//if ((LeadCount == 1 || LeadCount == 2) && ((LeadCount - *LeadBeg - *(LeadEnd - 1)) > 0))
		if (LeadCount > 0 && LeadCount < maxSigBins)
		{
			SigArray[2] = true;
			SigLeadPosition[1] = std::find(LeadBeg, LeadEnd, 1);
		}

		//if ((LagCount == 1 || LagCount == 2) && ((LagCount - *LagBeg - *(LagEnd - 1)) > 0))
		if (LagCount > 0 && LagCount < maxSigBins)
		{
			SigArray[3] = true;
			SigLagPosition[1] = std::find(std::make_reverse_iterator(LagEnd), std::make_reverse_iterator(LagBeg + 1), 1);
		}


		//Verify that Lead Excitatory Correlation and Lead Inhibitory Correlation exists.
		//Keeps the peak or trough that is closer to the reference spikes.
		if (SigArray[2] && SigArray[0])
		{
			if (SigLeadPosition[0] < SigLeadPosition[1])
			{
				SigArray[2] = false;
			}
			else
			{
				SigArray[0] = false;
			}
		}

		if (SigArray[3] && SigArray[1])
		{
			if (SigLagPosition[0] < SigLagPosition[1])
			{
				SigArray[3] = false;
			}
			else
			{
				SigArray[1] = false;
			}
		}


	}

private:

	Experiment OdorEx;
	BrainRegion Reference;
	BrainRegion Target;

	std::vector<std::vector<uint32_t>> StimLockedSpikesRef;
	std::vector<std::vector<uint32_t>> StimLockedSpikesTar;

	std::atomic<int> BinSize;
	std::atomic<int> Epoch;
	std::atomic<int> NoBins;
	std::atomic<int> GlobalReferenceUnit = 0;
	std::atomic<int> GlobalTargetUnit = 0;

	std::random_device Rd;
	std::default_random_engine Generator;
	//std::mt19937 Generator;

	std::mutex muVars;
	std::mutex muios;
	std::mutex muFile;
	//std::mutex muGenerator;

	std::ofstream CorrFile;
	//std::ofstream JitteredMatrixFile("JitteredMatrix" + std::to_string(Stimulus + 1) + ".txt");


};

constexpr unsigned char SHUFFLING = 0;
constexpr unsigned char BASICJITTER = 1;
constexpr unsigned char INTERJITTER = 2;

constexpr unsigned char ZTEST = 0;
constexpr unsigned char PERMUTATIONTEST = 1;

#define STALowerBoundTExc(FI, LI, Limit) STAUpperBoundT(FI, LI, Limit);
#define STAUpperBoundTExc(FI, LI, Limit) STALowerBoundT(FI, LI, Limit);