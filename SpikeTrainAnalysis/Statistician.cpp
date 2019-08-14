#include "Statistician.h"
#include <algorithm>
#include <future>
#include <iostream>
#include <numeric>
#include <cmath>
#include <string>
#include <utility>

#include <chrono>



Statistician::Statistician(std::string FileName, int BinSize, int Epoch, bool IsSpontaneous)
	:
	BinSize(BinSize),
	Epoch(Epoch),
	NoBins((Epoch / BinSize) * 2),
	BinSizeSec((double)BinSize / 1000.0),
	EpochSec((double)Epoch / 1000.0),
	OdorEx(FileName,IsSpontaneous),
	Reference(OdorEx.RDataFile(), OdorEx.GetUnitsRef(), OdorEx.GetRefSizePos(), OdorEx.GetRefTrainPos()),
	Target(OdorEx.RDataFile(), OdorEx.GetUnitsTar(), OdorEx.GetTarSizePos(), OdorEx.GetTarTrainPos()),
	StimLockedSpikesRef((long long)OdorEx.GetStimuli() * (long long)OdorEx.GetMagnitudes() * (long long)OdorEx.GetTrials() * (long long)OdorEx.GetUnitsRef()),
	StimLockedSpikesTar((long long)OdorEx.GetStimuli() * (long long)OdorEx.GetMagnitudes() * (long long)OdorEx.GetTrials() * (long long)OdorEx.GetUnitsTar()),
	Generator(Rd())
{
	if(IsSpontaneous)
	{ 
		StimLockedSpikesRef = Reference.UnitsCopy();
		StimLockedSpikesTar = Target.UnitsCopy();
	}
	else
	{
		SetStimLockedSpikes();
	}


}

Statistician::Statistician(std::string FileName, int BinSize, int Epoch, double Interval)
	:
	BinSize(BinSize),
	Epoch(Epoch),
	NoBins((Epoch / BinSize) * 2),
	BinSizeSec((double)BinSize / 1000.0),
	EpochSec((double)Epoch / 1000.0),
	OdorEx(FileName,false),
	Reference(OdorEx.RDataFile(), OdorEx.GetUnitsRef(), OdorEx.GetRefSizePos(), OdorEx.GetRefTrainPos()),
	Target(OdorEx.RDataFile(), OdorEx.GetUnitsTar(), OdorEx.GetTarSizePos(), OdorEx.GetTarTrainPos()),
	StimLockedSpikesRef((long long)OdorEx.GetStimuli() * (long long)OdorEx.GetMagnitudes() * (long long)OdorEx.GetTrials() * (long long)OdorEx.GetUnitsRef()),
	StimLockedSpikesTar((long long)OdorEx.GetStimuli() * (long long)OdorEx.GetMagnitudes() * (long long)OdorEx.GetTrials() * (long long)OdorEx.GetUnitsTar()),
	Generator(Rd())
{
	SetPREXLockedSpikes(Interval);
}

void Statistician::SetStimLockedSpikes()
{
	auto ItSLSpikesRef = StimLockedSpikesRef.begin(); //Iterator to the begining of the Vector.
	auto ItSLSpikesTar = StimLockedSpikesTar.begin(); //Iterator to the begining of the Vector.

	for (auto On = OdorEx.GetStimOn().begin(), Off = OdorEx.GetStimOff().begin(), end = OdorEx.GetStimOn().end();
		On < end;
		++On, ++Off)
	{
		for (auto Unit = Reference.RUnits().begin(), LastUnit = Reference.RUnits().end();
			Unit < LastUnit;
			++Unit)
		{
			ItSLSpikesRef->resize(Unit->size()); //Initializing the wraped vectors with the size of the total train

			//Using std algorithms library to get just the units in the interval [On Off].
			ItSLSpikesRef->erase(std::copy_if(Unit->begin(),
				Unit->end(),
				ItSLSpikesRef->begin(),
				[On, Off](double SpikeTime) {return SpikeTime > (*On) && SpikeTime < (*Off); }) //Lambda as predicate for the algorithm.
				,ItSLSpikesRef->end()); 

			ItSLSpikesRef->shrink_to_fit(); //Free garbage memory of the Vector. Very important cause its HUGE waste in this case.

			++ItSLSpikesRef;
		}

		for (auto Unit = Target.RUnits().begin(), LastUnit = Target.RUnits().end();
			Unit < LastUnit;
			++Unit)
		{
			ItSLSpikesTar->resize(Unit->size()); //Initializing the wraped vectors with the size of the total train

			//Using std algorithms library to get just the units in the interval [On Off].
			ItSLSpikesTar->erase(std::copy_if(Unit->begin(),
				Unit->end(),
				ItSLSpikesTar->begin(),
				[On, Off](double SpikeTime) {return SpikeTime > (*On) && SpikeTime < (*Off); }) //Lambda as predicate for the algorithm.
				, ItSLSpikesTar->end());

			ItSLSpikesTar->shrink_to_fit(); //Free garbage memory of the Vector. Very important cause its HUGE waste in this case.

			++ItSLSpikesTar;
		}
	}
}

void Statistician::SetPREXLockedSpikes(double Interval)
{
	auto ItSLSpikesRef = StimLockedSpikesRef.begin(); //Iterator to the begining of the Vector.
	auto ItSLSpikesTar = StimLockedSpikesTar.begin(); //Iterator to the begining of the Vector.

	for (auto PREXOn = OdorEx.GetPREXTimes().begin(), end = OdorEx.GetPREXTimes().end();
		PREXOn < end;
		++PREXOn)
	{
		for (auto Unit = Reference.RUnits().begin(), LastUnit = Reference.RUnits().end();
			Unit < LastUnit;
			++Unit)
		{
			ItSLSpikesRef->resize(Unit->size()); //Initializing the wraped vectors with the size of the total train

			//Using std algorithms library to get just the units in the interval [On Off].
			ItSLSpikesRef->erase(std::copy_if(Unit->begin(),
				Unit->end(),
				ItSLSpikesRef->begin(),
				[PREXOn,Interval](double SpikeTime) {return SpikeTime > (*PREXOn) && SpikeTime < (*PREXOn) + Interval; }) //Lambda as predicate for the algorithm.
				, ItSLSpikesRef->end());

			ItSLSpikesRef->shrink_to_fit(); //Free garbage memory of the Vector. Very important cause its HUGE waste in this case.

			++ItSLSpikesRef;
		}

		for (auto Unit = Target.RUnits().begin(), LastUnit = Target.RUnits().end();
			Unit < LastUnit;
			++Unit)
		{
			ItSLSpikesTar->resize(Unit->size()); //Initializing the wraped vectors with the size of the total train

			//Using std algorithms library to get just the units in the interval [On Off].
			ItSLSpikesTar->erase(std::copy_if(Unit->begin(),
				Unit->end(),
				ItSLSpikesTar->begin(),
				[PREXOn,Interval](double SpikeTime) {return SpikeTime > (*PREXOn) && SpikeTime < (*PREXOn) + Interval; }) //Lambda as predicate for the algorithm.
				, ItSLSpikesTar->end());

			ItSLSpikesTar->shrink_to_fit(); //Free garbage memory of the Vector. Very important cause its HUGE waste in this case.

			++ItSLSpikesTar;
		}
	}
}

void Statistician::SpikeTrainCorr(const std::vector<double>& reference, const std::vector<double>& target, std::vector<unsigned int>& Spikes, unsigned int& Count)
{
	//Computes Correlations separated by bins;

	double CurrentBinF;
	double CurrentBinL;
	auto LBit = target.begin();

		//std::cout << "Shuff: " << "\n";
	//for(auto Spike = reference.begin(), LastSpike = reference.end(); Spike < LastSpike; ++Spike)
	for (const double& Spike : reference)
	{
		CurrentBinF = Spike - EpochSec; // Set the current bins for the lambda function.
		CurrentBinL = CurrentBinF + BinSizeSec;

		//Boundaries of the target spikes.
		LBit = std::lower_bound(LBit, target.end(), Spike - EpochSec);
		auto UBit = std::upper_bound(LBit, target.end(), Spike + EpochSec);

		//Ierators for Count Corr vec
		auto Bin = Spikes.begin(), LastBin = Spikes.end();

		for (; Bin < LastBin - (NoBins/2); ++Bin)
		{
			*Bin += (unsigned int)std::count_if(LBit, //std library stuff, very convenient and fast.
				UBit,
				[&CurrentBinF, &CurrentBinL](double TargetSpike) 
				{
					return TargetSpike >= CurrentBinF && TargetSpike < CurrentBinL; 
				}
			);
			CurrentBinF = CurrentBinL;
			CurrentBinL = CurrentBinF + BinSizeSec;
		}

		for (; Bin < LastBin; ++Bin)
		{
			*Bin += (unsigned int)std::count_if(LBit, //std library stuff, very convenient and fast.
				UBit,
				[&CurrentBinF, &CurrentBinL](double TargetSpike)
				{
					return TargetSpike > CurrentBinF && TargetSpike <= CurrentBinL;
				}
			);
			CurrentBinF = CurrentBinL;
			CurrentBinL = CurrentBinF + BinSizeSec;
		}
	}
		//std::cout << "Shuff2: " << "\n";
	Count += (unsigned int)reference.size();
}

void Statistician::SpikeTrainJitter(const std::vector<double>& reference, const std::vector<double>& target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count)
{

	//This is bad because it is just copy pasta from the SpikeTrainCorr func and adding just a single Line for Jittering short intervals, but Im lazy to change the arquitecture.
	//I think that passing a bool to SpikeTrainCorr would just make the code messier.


	//Setting Pseudo Random Number Uniform Distribution for jittering [-5ms, 5ms] (Fujisawa, 2018)
	std::uniform_real_distribution<double> distribution(-0.005, 0.005);
	std::vector<double> JitteredTarget(target.size());


	for (auto Spikes = SpikesMatrix.begin(), SMEnd = SpikesMatrix.end(); Spikes < SMEnd; ++Spikes)
	{
		//Computes Correlations separated by bins;

		//std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
		std::transform(target.cbegin(), target.cend(), JitteredTarget.begin(), //Jittering here!
			[this, &distribution](const double& Spike) { return Spike + distribution(Generator); });


		SpikeTrainCorr(reference, JitteredTarget, *Spikes,Count);

		/*std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::chrono::duration<float> duration = end - start;
		std::cout << duration.count() << "\n";*/
	}

}

void Statistician::SpikeTrainShuffle(const std::vector<double>& reference, std::vector<double> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count)
{

	double TargetMin;
	double TargetMax;
	double TargetRandom;

	//Setting Pseudo Random Number Uniform Distribution
	std::uniform_int_distribution<int> distribution(0, (int)target.size()-1);

	for (auto Spikes = SpikesMatrix.begin(), SMEnd = SpikesMatrix.end(); Spikes < SMEnd ; ++Spikes)
	{
		//The original cross corr is included in the actual random distribution, jittering seems like a better method.
		auto RandomIt = target.begin() + distribution(Generator);
		TargetMin = *target.begin();
		TargetMax = *(target.end()-1);
		TargetRandom = *RandomIt;

		std::for_each(target.begin(), 
			RandomIt,
			[&TargetMax, &TargetMin, &TargetRandom](double& Spike) { Spike += (TargetMax + TargetMin - TargetRandom); });

		std::for_each(RandomIt,
			target.end(),
			[&TargetMin, &TargetRandom](double& Spike) { Spike -= (TargetRandom - TargetMin); });

		std::sort(target.begin(), target.end());

		SpikeTrainCorr(reference,
			target,
			*Spikes,
			Count);
	}
}

void Statistician::RunSingleThread(int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double ZorPVal, bool ExcZeroLag)
{
	for (int Stimulus = 0; Stimulus < OdorEx.GetStimuli(); Stimulus++)
	{
		MasterSpikeCrossCorrWorker(Stimulus, ResampledSets, ResamplingMethod, StatTest, ZorPVal, ExcZeroLag);
	}
}

void Statistician::RunThreadPool(int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double ZorPVal, bool ExcZeroLag)
{
	//ThreadPool with n threads. n = number of odors.
	std::vector<std::future<void>> ThreadPool(OdorEx.GetStimuli());

	auto CurrentThread = ThreadPool.begin();
	auto EndThread = ThreadPool.end();
	
	for (int Stimulus = 0; CurrentThread < EndThread; ++CurrentThread, Stimulus++)
	{
		*CurrentThread = std::async(std::launch::async, &Statistician::MasterSpikeCrossCorrWorker,
			this, Stimulus, ResampledSets, ResamplingMethod, StatTest, ZorPVal, ExcZeroLag);
	}

	while (true)
	{
		//This while loop freezes at the end, I need to implement in another way, but my thread pool works as the single thread.
		if (std::all_of(ThreadPool.begin(), ThreadPool.end(),
			[](std::future<void>& CurrentThread) {return CurrentThread.valid(); }))
		{
			for (CurrentThread = ThreadPool.begin(); CurrentThread < EndThread; ++CurrentThread)
			{
				CurrentThread->get();
			}
			return;
		}
	}
}

void Statistician::MasterSpikeCrossCorrWorker(int Stimulus, int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double PVal, bool ExcZeroLag)
{

	//NOTE: Im not convienced that STD matrices are the best way to deal with the problem. They are well allocated but anyway they may impact the performance,
	//STD problem may be solved with the use of other statistics instead of Z test (Fujisawa, 2008).
	//it is a posibility to implement fujisawa statistics but I need to try them first on MATLAB.

	mu.lock();

	//Locked code to access common memory between threads
	int UnitsRef = OdorEx.GetUnitsRef();
	int UnitsTar = OdorEx.GetUnitsTar();
	int Trials = OdorEx.GetTrials();
	auto SLSRB = StimLockedSpikesRef.cbegin();
	auto SLSTB = StimLockedSpikesTar.cbegin();

	std::cout << "Stimulus: " << Stimulus +1 << ", Ref: " << UnitsRef 
		<< ", Tar: " << UnitsTar << ", Trials: " << Trials << "\n";

	//Put this thread to sleep just for debugging puposes.
	//std::this_thread::sleep_for(std::chrono::milliseconds(1000));
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<unsigned int> SpikesCountCorr(NoBins); //Main raw correlation Vector

	std::vector<std::vector<unsigned int>> SpikesCountResampled(ResampledSets,std::vector<unsigned int>(NoBins)); // Good! Resampling Matrix, this is annoying but necessary to obtain the standard deviation.
	std::vector<unsigned int> SpikesSTDCount(ResampledSets);

	//Vars when working with PermTest comp fujisawa, 2008.
	std::vector<uint32_t> LPWBand(NoBins); //
	std::vector<uint32_t> UPWBand(NoBins);
	int PValPlace = (int)std::ceil(double(SpikesCountResampled.size() * PVal) / 2.0);
	std::vector<std::vector<uint32_t>> LPWBands(PValPlace, std::vector<uint32_t>(NoBins));
	std::vector<std::vector<uint32_t>> UPWBands(PValPlace, std::vector<uint32_t>(NoBins));
	std::pair<uint32_t, uint32_t> GlobalBands(0, 0);
	mu.unlock();

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//File for storing the sig data.
	//std::ofstream CorrFile("Stimulus" + std::to_string(Stimulus + 1) + ".SigCorr", std::ios::binary);
	std::ofstream CorrFile("Stimulus" + std::to_string(Stimulus + 1) + ".txt");


	//Check if we want to exclude "zero lag" correlations.
	int BinExcluded = 0;
	if (ExcZeroLag)
		BinExcluded = 1;
	//Nested loops for running the whole analysis. There may be some improvement specially in the las loop if the data is parsed better from matlab.

	//Stimulus locked reference spike train loop
	unsigned short ReferenceUnit = 1;
	for (auto RefTrain = SLSRB + ((__int64)Stimulus * UnitsRef),
		endRT = RefTrain + UnitsRef;
		RefTrain < endRT
		; ++RefTrain, ReferenceUnit++)
	{
		if (ReferenceUnit == 1)
		{
			//Stimulus locked target spike train loop
			unsigned short TargetUnit = 1;
			for (auto TarTrain = SLSTB + ((__int64)Stimulus * UnitsTar),
				endTT = TarTrain + UnitsTar;
				TarTrain < endTT
				; ++TarTrain, TargetUnit++)
			{
				auto RefTrialTrain = RefTrain; //this is the downside of the way I parse the matlab data.
				auto TarTrialTrain = TarTrain; //Aux vars to prevent modification of original vars.
				unsigned int CountCorr = 0;
				unsigned int CountRes = 0;

				//Trial Loop/////////////////////////////////////////////////////////////////////////
				for (int Trial = 0; Trial < Trials;
					RefTrialTrain += UnitsRef, TarTrialTrain += UnitsTar, Trial++)
				{
					if ((RefTrialTrain->size() != 0 && TarTrialTrain->size() != 0)) //Check if trains are not empty.
					{
						SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, CountCorr); //Compute Corr.

						switch (ResamplingMethod)
						{
						case SHUFFLING:
							SpikeTrainShuffle(*RefTrialTrain, *TarTrialTrain, SpikesCountResampled, CountRes); //Compute Corr shuffling method.
							break;

						case JITTERING:
							SpikeTrainJitter(*RefTrialTrain, *TarTrialTrain, SpikesCountResampled, CountRes); //Compute Corr Jittering method.
							break;

						default:
							SpikeTrainShuffle(*RefTrialTrain, *TarTrialTrain, SpikesCountResampled, CountRes);
							break;
						}
					}
				}
				/////////////////////////////////////////////////////////////////////////////////////

				//////////////////////////////////////////Use to be PrepPermTest/////////////////////////////////////////////////////////////////////////////////


				//Mean and STD of the matrix and vectors of the choosen resampling method.///////////
				//CountRes /= ResampledSets; //This needs to be divided into ResampledSets because that is the size of the Matrix, is not a vector anymore.
				bool GoodData = true;
				bool GoodAlpha = false;

				for (int Bin = 0; Bin < NoBins; Bin++)
				{
					//Looping through the Matrix and filling the STDCount Vector.
					auto STDCount = SpikesSTDCount.begin();
					auto STDCountEnd = SpikesSTDCount.end();
					for (auto BinVec = SpikesCountResampled.cbegin(), BinVecEnd = SpikesCountResampled.cend();
						BinVec < BinVecEnd;
						++BinVec, ++STDCount)
					{
						*STDCount = *(BinVec->begin() + Bin);

					}

					//If we have a zero this means we have an unpopulated resampled correlogram which is useless for statistical comparisons.
					if (std::accumulate(SpikesSTDCount.begin(), SpikesSTDCount.end(), 0) == 0)
					{
						GoodData = false;
						break;
					}

					//Sorting the Resampled data to get the points at the desire PVal
					std::sort(SpikesSTDCount.begin(), SpikesSTDCount.end());

					//Filling the Pointwise bands Matrix.
					int ProvPlace = PValPlace;
					for (auto LPWsit = LPWBands.begin(), UPWsit = UPWBands.begin(), End = LPWBands.end();
						LPWsit < End; ++LPWsit, ++UPWsit)
					{

						*(LPWsit->begin() + Bin) = *(SpikesSTDCount.begin() + ProvPlace - 1);
						*(UPWsit->begin() + Bin) = *(SpikesSTDCount.end() - ProvPlace);

						ProvPlace--;
					}

					//Pointwise bands.
					*(LPWBand.begin() + Bin) = *(SpikesSTDCount.begin() + PValPlace - 1); // Low Pval.
					*(UPWBand.begin() + Bin) = *(SpikesSTDCount.end() - PValPlace); // Upper Val.

				}

				//Loop for gettting the P of surrogate data sets that break the Pointwise bands at ANY point, that is the PVal of the global band that corresponds to the alpha of the Pairwise bands.
				for (auto LPWVecit = LPWBands.cbegin(), UPWVecit = UPWBands.cbegin(), Ends = LPWBands.cend();
					LPWVecit < Ends; ++LPWVecit, ++UPWVecit)
				{
					uint32_t PWCount = 0;

					for (auto BinVec = SpikesCountResampled.cbegin(), BinVecEnd = SpikesCountResampled.cend();
						BinVec < BinVecEnd;
						++BinVec)
					{
						for (auto LPWit = LPWVecit->cbegin(), UPWit = UPWVecit->cbegin(), End = LPWVecit->cend(), ResDatait = BinVec->cbegin();
							LPWit < End; ++LPWit, ++UPWit, ++ResDatait)
						{
							if (*ResDatait < *LPWit || *ResDatait > * UPWit)
							{
								PWCount += 1;
								break;
							}
						}
					}

					if ((double)PWCount / (double)SpikesCountResampled.size() <= PVal / 2.0)
					{
						//Defining global bands.
						auto LowBand = std::min_element(LPWVecit->begin(), LPWVecit->end());
						auto UpperBand = std::max_element(UPWVecit->begin(), UPWVecit->end());

						GlobalBands.first = *LowBand;
						GlobalBands.second = *UpperBand;

						GoodAlpha = true;
					}
				}


				if (!GoodAlpha)
				{
					//Defining global bands as the extremes in case the Pval is too low.
					auto LowBand = std::min_element(LPWBands.rbegin()->begin(), LPWBands.rbegin()->end());
					auto UpperBand = std::max_element(UPWBands.rbegin()->begin(), UPWBands.rbegin()->end());

					GlobalBands.first = *LowBand;
					GlobalBands.second = *UpperBand;

				}

				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				/////////////////////////////////////////////Permutation Test///////////////////////////////////////////////////////////////////////////////////

				if ((CountCorr != 0) && GoodData)
				{
					//Writing Significant Correlations to File.

					bool LeadEx = false;
					bool LagEx = false;
					bool LeadIn = false;
					bool LagIn = false;

					//Check significant corrs.
					auto LeadExC = std::count_if(SpikesCountCorr.end() - (SpikesCountCorr.size() / 2), SpikesCountCorr.end(),
						[&GlobalBands](uint32_t& RawVal)
						{
							return RawVal > GlobalBands.second;
						});
					auto LagExC = std::count_if(SpikesCountCorr.begin(), SpikesCountCorr.begin() + (SpikesCountCorr.size() / 2),
						[&GlobalBands](uint32_t& RawVal)
						{
							return RawVal > GlobalBands.second;
						});
					auto LeadInC = std::count_if(SpikesCountCorr.end() - (SpikesCountCorr.size() / 2), SpikesCountCorr.end(),
						[&GlobalBands](uint32_t& RawVal)
						{
							return RawVal < GlobalBands.first;
						});
					auto LagInC = std::count_if(SpikesCountCorr.begin(), SpikesCountCorr.begin() + (SpikesCountCorr.size() / 2),
						[&GlobalBands](uint32_t& RawVal)
						{
							return RawVal < GlobalBands.first;
						});

					//Verify that sig corrs are not due to common stimulus modulation. thin significant bins represent that.
					if (LeadExC == 1 || LeadExC == 2)
						LeadEx = true;

					if (LagExC == 1 || LagExC == 2)
						LagEx = true;

					if (LeadInC == 1 || LeadInC == 2)
						LeadIn = true;

					if (LagInC == 1 || LagInC == 2)
						LagIn = true;

					mu.lock();
					if (LeadEx && LagEx)
					{
						//Code to store in txt files.
						CorrFile << 1 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}
					else if (LeadEx)
					{
						CorrFile << 2 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}
					else if (LagEx)
					{
						CorrFile << 3 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}

					if (LeadIn && LagIn)
					{
						CorrFile << 4 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}
					else if (LeadIn)
					{
						CorrFile << 5 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}
					else if (LagIn)
					{
						CorrFile << 6 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}

					if ((LeadIn || LagIn) && (LeadEx || LagEx))
					{
						CorrFile << 7 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}
					mu.unlock();
				}

				////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

				//Reseting Count Vectors and Matrix;
				std::fill(SpikesCountCorr.begin(), SpikesCountCorr.end(), 0);

				std::for_each(SpikesCountResampled.begin(), SpikesCountResampled.end(),
					[](std::vector<unsigned int>& BinVec)
					{
						std::fill(BinVec.begin(), BinVec.end(), 0);
					});

				mu.lock();
				std::cout << "Stimulus " << Stimulus + 1 << ". Finished reference unit " << ReferenceUnit << " vs target unit " << TargetUnit << ".\n";
				mu.unlock();
			}
		}
	}

	if (CorrFile.bad())
		std::cout << "bad";

	else if (CorrFile.eof())
		std::cout << "eof";

	else if (CorrFile.fail())
		std::cout << "other fail";

	else if (CorrFile.good())
	{
		CorrFile.close();
		mu.lock();
		std::cout << "Data file was closed successfully\n";
		mu.unlock();
	}

	if (CorrFile.rdstate() == (std::ios_base::failbit | std::ios_base::eofbit))
	{
		mu.lock();
		std::cout << "stream state is eofbit\n";
		mu.unlock();
	}
}