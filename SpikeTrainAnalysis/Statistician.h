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
	void SpikeTrainJitter(const std::vector<double>& reference, std::vector<double> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count);
	void SpikeTrainShuffle(const std::vector<double>& reference, std::vector<double> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count);
	void MasterSpikeCrossCorrWorker(int Stimulus, int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double ZorPVal, bool ExcZeroLag);
	void WriteToFileWorker(std::ofstream& CorrFile, std::vector<double>& CorrVec, uint32_t CorrCount);
	bool PrepZTest(std::vector<unsigned int>& SpikesSTDCount, std::vector<std::vector<unsigned int>>& SpikesCountResampled, std::vector<double>& SpikesSTDResampled, std::vector<double>& SpikesPResampled, uint32_t CountRes);
	bool PrepPermTest(std::vector<uint32_t>& SpikesSTDCount, std::vector<std::vector<uint32_t>>& SpikesCountResampled, std::vector<uint32_t>& LPWBand, std::vector<uint32_t>& UPWBand, std::pair<uint32_t, uint32_t>& GlobalBands, double PVal, int PValPlace, std::vector<std::vector<uint32_t>>& LPWBands, std::vector<std::vector<uint32_t>>& UPWBands);
	void ZTestToFile(std::vector<double>& SpikesSTDResampled, std::vector<unsigned int>& SpikesCountCorr, std::vector<double>& SpikesPCorr, std::vector<double>& SpikesPResampled, unsigned int CountCorr, int BinExcluded, double ZorPVal, std::ofstream& CorrFile, uint16_t ReferenceUnit, uint16_t TargetUnit, uint32_t CountRes);
	void PermTestToFile(std::pair<uint32_t, uint32_t>& GlobalBands, std::vector<uint32_t>& SpikesCountCorr, std::vector<uint32_t>& LPWBand, std::vector<uint32_t>& UPWBand, unsigned int CountCorr, std::ofstream& CorrFile, uint16_t ReferenceUnit, uint16_t TargetUnit);

	template <typename T>
	void WriteToFileWorkerT(std::ofstream& CorrFile, std::vector<T>& CorrVec, uint32_t CountCorr)
	{
		std::for_each(CorrVec.begin(), CorrVec.end(),
			[&CorrFile, &CountCorr](T& Bin)
			{
				CorrFile << Bin * CountCorr << ", ";
			});
		CorrFile << CountCorr << ", ";
	}

	//void SpikeTrainShift(); // I dont know if Im gonna implement shift, seems that is not very useful for my actual experiment.

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