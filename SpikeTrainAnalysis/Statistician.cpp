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
	OdorEx(FileName,IsSpontaneous),
	Reference(OdorEx.RDataFile(), OdorEx.GetUnitsRef(), OdorEx.GetRefSizePos(), OdorEx.GetRefTrainPos()),
	Target(OdorEx.RDataFile(), OdorEx.GetUnitsTar(), OdorEx.GetTarSizePos(), OdorEx.GetTarTrainPos()),
	StimLockedSpikesRef((long long)OdorEx.GetStimuli() * (long long)OdorEx.GetMagnitudes() * (long long)OdorEx.GetTrials() * (long long)OdorEx.GetUnitsRef()),
	StimLockedSpikesTar((long long)OdorEx.GetStimuli() * (long long)OdorEx.GetMagnitudes() * (long long)OdorEx.GetTrials() * (long long)OdorEx.GetUnitsTar()),
	Generator(Rd()),
	CorrFile("CorrFile.txt")

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

Statistician::Statistician(std::string FileName, int BinSize, int Epoch, uint32_t Interval)
	:
	BinSize(BinSize),
	Epoch(Epoch),
	NoBins((Epoch / BinSize) * 2),
	OdorEx(FileName,false),
	Reference(OdorEx.RDataFile(), OdorEx.GetUnitsRef(), OdorEx.GetRefSizePos(), OdorEx.GetRefTrainPos()),
	Target(OdorEx.RDataFile(), OdorEx.GetUnitsTar(), OdorEx.GetTarSizePos(), OdorEx.GetTarTrainPos()),
	StimLockedSpikesRef((long long)OdorEx.GetStimuli() * (long long)OdorEx.GetMagnitudes() * (long long)OdorEx.GetTrials() * (long long)OdorEx.GetUnitsRef()),
	StimLockedSpikesTar((long long)OdorEx.GetStimuli() * (long long)OdorEx.GetMagnitudes() * (long long)OdorEx.GetTrials() * (long long)OdorEx.GetUnitsTar()),
	Generator(Rd()),
	CorrFile("CorrFile.txt")
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
				[On, Off](uint32_t SpikeTime) {return SpikeTime > (*On) && SpikeTime < (*Off); }) //Lambda as predicate for the algorithm.
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
				[On, Off](uint32_t SpikeTime) {return SpikeTime > (*On) && SpikeTime < (*Off); }) //Lambda as predicate for the algorithm.
				, ItSLSpikesTar->end());

			ItSLSpikesTar->shrink_to_fit(); //Free garbage memory of the Vector. Very important cause its HUGE waste in this case.

			++ItSLSpikesTar;
		}
	}
}

void Statistician::SetPREXLockedSpikes(uint32_t Interval)
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
				[PREXOn,Interval](uint32_t SpikeTime) {return SpikeTime > (*PREXOn) && SpikeTime < (*PREXOn) + Interval; }) //Lambda as predicate for the algorithm.
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
				[PREXOn,Interval](uint32_t SpikeTime) {return SpikeTime > (*PREXOn) && SpikeTime < (*PREXOn) + Interval; }) //Lambda as predicate for the algorithm.
				, ItSLSpikesTar->end());

			ItSLSpikesTar->shrink_to_fit(); //Free garbage memory of the Vector. Very important cause its HUGE waste in this case.

			++ItSLSpikesTar;
		}
	}
}

void Statistician::SpikeTrainCorr(const std::vector<uint32_t>& reference, const std::vector<uint32_t>& target, std::vector<unsigned int>& Spikes, unsigned int& Count)
{
	//Computes Correlations separated by bins;

	long long CurrentBinF;
	long long CurrentBinL;
	auto LBit = target.begin();
	auto UBit = target.begin();

	//std::cout << "Shuff: " << "\n";
    //for(auto Spike = reference.begin(), LastSpike = reference.end(); Spike < LastSpike; ++Spike)
	for (const uint32_t& Spike : reference)
	{
		CurrentBinF = (long long)Spike - Epoch; // Set the current bins for the lambda function.
		CurrentBinL = CurrentBinF + BinSize;

		//Boundaries of the target spikes.

		STALowerBoundT(LBit, target.end(), CurrentBinF);
		STAUpperBoundT(UBit, target.end(), Spike + Epoch);

		auto First = LBit;
		auto Last = LBit;
		STAUpperBoundTExc(Last, UBit, CurrentBinL);

		//Ierators for Count Corr vec
		auto Bin = Spikes.begin(), LastBin = Spikes.end();



		/////This loops are written this way to avoid counting zero lag correlations. They are implemented using pointer
		// arithmetic with custom made functions. STA stands for Spike Train Analysis.

		for (; Bin < LastBin - (NoBins / 2) - 1 ; ++Bin)
		{
			*Bin += (unsigned int)std::distance(First, Last);
			CurrentBinF = CurrentBinL;
			CurrentBinL = CurrentBinF + BinSize;

			First = Last;
			STAUpperBoundTExc(Last, UBit, CurrentBinL);
		}

		*Bin += (unsigned int)std::distance(First, Last);
		CurrentBinF = CurrentBinL;
		CurrentBinL = CurrentBinF + BinSize;

		STALowerBoundTExc(First, UBit, CurrentBinF);
		STAUpperBoundT(Last, UBit, CurrentBinL);

		++Bin;

		for (; Bin < LastBin; ++Bin)
		{
			*Bin += (unsigned int)std::distance(First, Last);
			CurrentBinF = CurrentBinL;
			CurrentBinL = CurrentBinF + BinSize;

			First = Last;
			STAUpperBoundT(Last, UBit, CurrentBinL);
		}
	}
	Count += (unsigned int)reference.size();
}

void Statistician::SpikeTrainIntervalJitter(const std::vector<unsigned int>& Spikes, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count)
{

	// This is the correct and possibly final implementation for the interval jittering method.
	// Interval jittering provides good resampled data to test my H0.

	// H0 = All placements of spikes under JitterInterval are equally probable.
	// That is, there is no temporal structure at a scale finer than JitterInterval.
	// Remember that in order to test this H0 is very important to have the same spike counts in the resampled histograms (this func does that),
	// if we dont have the same counts, testing H0 is impossible and makes no sense.

	// If JitterInterval is 150 that finer scale is 5 ms. If JitterInterval is 90, scale is 3 ms.
	// That is because we are using samples not time stamps and our Sampling Frequency is 30kHz.
	// I choose 3 ms for now because I notice some troughs and co-modulation in spikes that are preserved under 5 ms, 
	// is common input possible at such finer scales? or is odor representation in PfCx that fast?

	// References for interval jittering: Amarasingham(2011), Amarasingham(2012) and Platkiewicz(2017).

	// NOTE: Interval jittering is a great resampling method for spontaneous activity, when measuring time-locked connectivity it should be used along with Trial shuffling.
	//If the data breaks the sig bands of both tests we can assume a monosynaptic interaction that depends on the stimulus.
	std::default_random_engine Generator(Rd());
	uint32_t JitterInterval = 15 * 30; //This is 90 because I want an interval of 3 ms. 150 = 5 ms.


	//uint32_t JitterCounts = (Epoch / JitterInterval) * 2;
	uint32_t JitterCounts = (Epoch * 2 ) / JitterInterval;
	std::vector<uint32_t> JitterIntC(JitterCounts);
	auto JICIt = JitterIntC.begin();
	const auto JICEnd = JitterIntC.end();
	std::vector<std::uniform_int_distribution<uint32_t>> Distributions(JitterCounts);
	auto DistsIt = Distributions.begin();

	uint32_t IntervalConstant = (JitterInterval / BinSize);

	uint32_t DistF = 0;
	uint32_t DistL = DistF + JitterInterval;

	auto SpikeCount = Spikes.begin();

	for (; JICIt < (JICEnd - JitterCounts/2); SpikeCount += IntervalConstant, ++JICIt, ++DistsIt)
	{
		*JICIt = std::accumulate(SpikeCount, SpikeCount + IntervalConstant, 0);
		*DistsIt = std::uniform_int_distribution<uint32_t>(DistF, DistL - 1);

		DistF = DistL;
		DistL = DistF + JitterInterval;
	}

	for (; JICIt < JICEnd; SpikeCount += IntervalConstant, ++JICIt, ++DistsIt)
	{
		*JICIt = std::accumulate(SpikeCount, SpikeCount + IntervalConstant, 0);
		*DistsIt = std::uniform_int_distribution<uint32_t>(DistF + 1, DistL);

		DistF = DistL;
		DistL = DistF + JitterInterval;
	}

	//Fake RefSample is any element of the set {x:x>Epoch};
	uint32_t FakeRefSample = Epoch;
	std::vector<uint32_t> FakeReference(1, FakeRefSample);
	std::vector<uint32_t> FakeTarget(std::accumulate(JitterIntC.begin(), JitterIntC.end(), 0));

	for (auto JittVec = SpikesMatrix.begin(), SMEnd = SpikesMatrix.end(); JittVec < SMEnd; ++JittVec)
	{
		JICIt = JitterIntC.begin();
		DistsIt = Distributions.begin();
		auto FTIt = FakeTarget.begin();

		for (; JICIt < JICEnd; ++JICIt, ++DistsIt)
		{
			auto FTEnd = FTIt + *JICIt;
			for (; FTIt < FTEnd; ++FTIt)
			{
				//muGenerator.lock();
				*FTIt = (*DistsIt)(Generator);
				//muGenerator.unlock();
			}
		}
		std::sort(FakeTarget.begin(), FakeTarget.end());
		SpikeTrainCorr(FakeReference, FakeTarget, *JittVec, Count);
	}
}

void Statistician::SpikeTrainBasicJitter(const std::vector<uint32_t>& reference, const std::vector<uint32_t>& target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count)
{
	//Basic jitter or Center jitter. this is just a heuristic, doesnt work to test the actual null hypothesis.
	//Do not use for statistics. (Amarasingham, 2011)

	//Setting Pseudo Random Number Uniform Distribution for jittering [-5ms, 5ms] (Fujisawa, 2018)
	//Fujisawa window is not useful for me, I chose a smalles time window [-2ms, 2ms].
	
	//150 samples correspond to 5ms at 30000 kHz.
	std::uniform_int_distribution<int> distribution(-150, 150);
	//std::normal_distribution<double> distribution(0,0.001);
	
	std::vector<uint32_t> JitteredTarget(target.size());
	std::ofstream DistFile("Dist.txt");

	for (auto Spikes = SpikesMatrix.begin(), SMEnd = SpikesMatrix.end(); Spikes < SMEnd; ++Spikes)
	{
		//Computes Correlations separated by bins

		std::transform(target.cbegin(), target.cend(), JitteredTarget.begin(), //Jittering here!
			[this, &distribution](const uint32_t& Spike) { return Spike + distribution(Generator); });

		//std::sort(JitteredTarget.begin(), JitteredTarget.end());

		SpikeTrainCorr(reference, JitteredTarget, *Spikes,Count);
		//WriteToFileWorkerT(DistFile, *Spikes);
		//DistFile << "\n";
	}

}

void Statistician::SpikeTrainBasicCuJitter(const std::vector<uint32_t>& reference, std::vector<uint32_t> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count)
{
	//Basic cumulative jitter or Center cumulative jitter. this is just a heuristic, doesnt work to test the actual null hypothesis.
	//Do not use for statistics.
	//Takes cumultive sums of random spike places.
	

	//Setting Pseudo Random Number Uniform Distribution for jittering [-5ms, 5ms] (Fujisawa, 2018)
	//Fujisawa window is not useful for me, I chose a smalles time window [-2ms, 2ms].
	
	//150 samples correspond to 5ms at 30000 kHz.
	std::uniform_int_distribution<int> distribution(-150, 150);
	//std::normal_distribution<double> distribution(0,0.001);

	std::ofstream DistFile("Dist.txt");

	for (auto Spikes = SpikesMatrix.begin(), SMEnd = SpikesMatrix.end(); Spikes < SMEnd; ++Spikes)
	{
		//Computes Correlations separated by bins

		std::for_each(target.begin(), target.end(),
			[this, &distribution](uint32_t& Spike) { Spike += distribution(Generator); });

		//std::sort(target.begin(), target.end());

		SpikeTrainCorr(reference, target, *Spikes, Count);

		//WriteToFileWorkerT(DistFile, *Spikes);
		//DistFile << "\n";
	}

}

void Statistician::SpikeTrainShuffle(const std::vector<uint32_t>& reference, std::vector<uint32_t> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count)
{

	uint32_t TargetMin;
	uint32_t TargetMax;
	uint32_t TargetRandom;

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
			[&TargetMax, &TargetMin, &TargetRandom](uint32_t& Spike) { Spike += (TargetMax + TargetMin - TargetRandom); });

		std::for_each(RandomIt,
			target.end(),
			[&TargetMin, &TargetRandom](uint32_t& Spike) { Spike -= (TargetRandom - TargetMin); });

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
		MasterSpikeCrossCorrSingleThread(Stimulus, ResampledSets, ResamplingMethod, StatTest, ZorPVal, ExcZeroLag);
	}
	CloseFiles();
}

void Statistician::RunThreadPool(int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double ZorPVal, bool ExcZeroLag)
{
	//ThreadPool with n threads. n = number of odors.

	unsigned int NoThreads = std::thread::hardware_concurrency();

	if (NoThreads < 4)
	{
		std::cout << "Multithreading not supported by the system "<< "\n";

	}
	else if (NoThreads > 16)
	{

		NoThreads = 16;

	}

	NoThreads = NoThreads - 2;

	std::vector<std::future<void>> ThreadPool(NoThreads);

	auto CurrentThread = ThreadPool.begin();
	auto EndThread = ThreadPool.end();
	
	for (int ThreadNo = 0; CurrentThread < EndThread; ++CurrentThread, ThreadNo++)
	{
		*CurrentThread = std::async(std::launch::async, &Statistician::MasterSpikeCrossCorrWorker,
			this, ThreadNo, ResampledSets, ResamplingMethod, StatTest, ZorPVal, ExcZeroLag);
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

			CloseFiles();

			return;
		}
	}
}

void Statistician::MasterSpikeCrossCorrWorker(int ThreadNo, int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double PVal, bool ExcZeroLag)
{

	//NOTE: Im not convienced that STD matrices are the best way to deal with the problem. They are well allocated but anyway they may impact the performance,
	//STD problem may be solved with the use of other statistics instead of Z test (Fujisawa, 2008).
	//it is a posibility to implement fujisawa statistics but I need to try them first on MATLAB.

	muVars.lock();

	

	//Locked code to access common memory between threads
	int UnitsRef = OdorEx.GetUnitsRef();
	int UnitsTar = OdorEx.GetUnitsTar();
	int Trials = OdorEx.GetTrials();


	std::cout << "ThreadNo: " << ThreadNo << ", Ref: " << UnitsRef
		<< ", Tar: " << UnitsTar << ", Trials: " << Trials << " ..... ";

	std::vector<std::vector<uint32_t>> SLSR4Thread(StimLockedSpikesRef);
	std::vector<std::vector<uint32_t>> SLST4Thread(StimLockedSpikesTar);

	auto SLSRB = SLSR4Thread.cbegin();
	auto SLSTB = SLST4Thread.cbegin();
	
	//Put this thread to sleep just for debugging puposes.
	//std::this_thread::sleep_for(std::chrono::milliseconds(1000));
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<unsigned int> SpikesCountCorr(NoBins); //Main raw correlation Vector

	std::vector<std::vector<unsigned int>> SpikesCountResampled(ResampledSets,std::vector<unsigned int>(NoBins)); // Good! Resampling Matrix, this is annoying but necessary to obtain the standard deviation.
	std::vector<unsigned int> SpikesSTDCount(ResampledSets);
	std::vector<unsigned int>AllSurrogateSpikeCounts(ResampledSets * NoBins);

	//Vars when working with PermTest comp fujisawa, 2008.
	std::vector<uint32_t> LPWBand(NoBins); //
	std::vector<uint32_t> UPWBand(NoBins);
	int PValPlace = (int)std::ceil(double(ResampledSets * PVal) / 2.0);
	std::vector<std::vector<uint32_t>> LPWBands(PValPlace, std::vector<uint32_t>(NoBins));
	std::vector<std::vector<uint32_t>> UPWBands(PValPlace, std::vector<uint32_t>(NoBins));
	std::pair<uint32_t, uint32_t> GlobalBands(0, 0);

	std::cout << " Thread Ready " << "\n";

	muVars.unlock();

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//File for storing the sig data.
	//std::ofstream CorrFile("Stimulus" + std::to_string(Stimulus + 1) + ".txt");

	//Check if we want to exclude "zero lag" correlations.
	int BinExcluded = 0;
	if (ExcZeroLag)
		BinExcluded = 1;
	
	//This apparently arbitrary numbers are the GlobalReferenceUnit and GlobalTargetUnit values before any modification. They
	//should be the same as in Statistician.h
	int PrevRefUnit = 0;
	int PrevTarUnit = 0;
	int ReferenceUnit = 0;
	int TargetUnit = 0;

	auto RefTrain = SLSRB; //+ ((__int64)Stimulus * UnitsRef); 

	auto TarTrain = SLSTB; //+ ((__int64)Stimulus * UnitsTar);

	//This loop is only inteded for Correlations whitin the same shank!!!!
	while (true)
	{
		muVars.lock();
		
		if (GlobalReferenceUnit == (UnitsRef - 2)) //&& GlobalTargetUnit == (UnitsTar - 1))
		{
			muVars.unlock();
			break;
		}
		else if (GlobalTargetUnit == (UnitsTar - 1))
		{
			GlobalReferenceUnit++;
			GlobalTargetUnit = GlobalReferenceUnit + 1;

			ReferenceUnit = GlobalReferenceUnit;
			TargetUnit = GlobalTargetUnit;


			RefTrain += ((long long)ReferenceUnit - (long long)PrevRefUnit);
				
			TarTrain = SLSTB + GlobalTargetUnit;// ((__int64)Stimulus * UnitsTar)


			PrevRefUnit = ReferenceUnit;
			PrevTarUnit = TargetUnit;

		}
		else
		{
			GlobalTargetUnit++;

			ReferenceUnit = GlobalReferenceUnit;
			TargetUnit = GlobalTargetUnit;

			RefTrain += ((long long)ReferenceUnit - (long long)PrevRefUnit);
			TarTrain += ((long long)TargetUnit - (long long)PrevTarUnit);

			PrevRefUnit = ReferenceUnit;
			PrevTarUnit = TargetUnit;
		}
	
		muVars.unlock();

		//if (TargetUnit == 2)
		if(ReferenceUnit < TargetUnit)
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
					switch (ResamplingMethod)
					{
					case SHUFFLING:
						SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, CountCorr); //Compute Corr.
						SpikeTrainShuffle(*RefTrialTrain, *TarTrialTrain, SpikesCountResampled, CountRes); //Compute Corr shuffling method.
						break;

					case BASICJITTER:
						SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, CountCorr);
						SpikeTrainBasicJitter(*RefTrialTrain, *TarTrialTrain, SpikesCountResampled, CountRes); //Compute Basic Jittering method.
						break;

					case INTERJITTER:
						SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, CountCorr);
						SpikeTrainIntervalJitter(SpikesCountCorr, SpikesCountResampled,CountRes);
						//SpikeTrainIntervalJitter4(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, SpikesCountResampled, CountCorr);
						//WriteToFileWorkerT(CorrFile, SpikesCountCorr); CorrFile << "\n";

						break;

					default:
						SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, CountCorr);
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
			auto AllCountIt = AllSurrogateSpikeCounts.begin();
			for (int Bin = 0; Bin < NoBins; Bin++)
			{
				//Looping through the Matrix and filling the STDCount Vector.
				auto STDCount = SpikesSTDCount.begin();
				auto STDCountEnd = SpikesSTDCount.end();
				for (auto BinVec = SpikesCountResampled.cbegin(), BinVecEnd = SpikesCountResampled.cend();
					BinVec < BinVecEnd;
					++BinVec, ++STDCount, ++AllCountIt)
				{
					*STDCount = *(BinVec->begin() + Bin);
					*AllCountIt = *(BinVec->begin() + Bin);

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
				/*int ProvPlace = PValPlace;
				for (auto LPWsit = LPWBands.begin(), UPWsit = UPWBands.begin(), End = LPWBands.end();
					LPWsit < End; ++LPWsit, ++UPWsit)
				{

					*(LPWsit->begin() + Bin) = *(SpikesSTDCount.begin() + ProvPlace - 1);
					*(UPWsit->begin() + Bin) = *(SpikesSTDCount.end() - ProvPlace);

					ProvPlace--;
				}*/

				//Pointwise bands.
				*(LPWBand.begin() + Bin) = *(SpikesSTDCount.begin() + PValPlace - 1); // Low Pval.
				*(UPWBand.begin() + Bin) = *(SpikesSTDCount.end() - PValPlace); // Upper Val.

			}

			//Sortiing ALL spike counts and getting the global bands.
			std::sort(AllSurrogateSpikeCounts.begin(), AllSurrogateSpikeCounts.end());

			GlobalBands.first = *(AllSurrogateSpikeCounts.begin() + PValPlace * 1 - 1);
			GlobalBands.second = *(AllSurrogateSpikeCounts.end() - PValPlace * 1);
			GoodAlpha = true;


			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			/////////////////////////////////////////////Permutation Test///////////////////////////////////////////////////////////////////////////////////

			if ((CountCorr != 0) && GoodData)
			{
				//Writing Significant Correlations to File.

				bool SigArray[4] = { false };

				GetSignifcantCorr(SpikesCountCorr, SigArray, GlobalBands, LPWBand, UPWBand);

				muFile.lock();
						
				if ((SigArray[2] || SigArray[3]) && (SigArray[0] || SigArray[1]))
				{
					CorrFile << 7 << ", " << ReferenceUnit + 1 << ", " << TargetUnit + 1 << ", ";
					WriteToFileWorkerT(CorrFile, SpikesCountCorr);
					WriteToFileWorkerT(CorrFile, LPWBand);
					WriteToFileWorkerT(CorrFile, UPWBand);
					CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
				}
				else
				{
					if (SigArray[0] && SigArray[1])
					{
						//Code to store in txt files.
						CorrFile << 1 << ", " << ReferenceUnit + 1 << ", " << TargetUnit + 1 << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}
					else if (SigArray[0])
					{
						CorrFile << 2 << ", " << ReferenceUnit + 1 << ", " << TargetUnit + 1 << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}
					else if (SigArray[1])
					{
						CorrFile << 3 << ", " << ReferenceUnit + 1 << ", " << TargetUnit + 1 << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}

					if (SigArray[2] && SigArray[3])
					{
						CorrFile << 4 << ", " << ReferenceUnit + 1 << ", " << TargetUnit + 1 << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}
					else if (SigArray[2])
					{
						CorrFile << 5 << ", " << ReferenceUnit + 1 << ", " << TargetUnit + 1 << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}
					else if (SigArray[3])
					{
						CorrFile << 6 << ", " << ReferenceUnit + 1 << ", " << TargetUnit + 1 << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
					}
				}


				muFile.unlock();
			}

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			//Reseting Count Vectors and Matrix;
			std::fill(SpikesCountCorr.begin(), SpikesCountCorr.end(), 0);

			std::for_each(SpikesCountResampled.begin(), SpikesCountResampled.end(),
				[](std::vector<unsigned int>& BinVec)
				{
					std::fill(BinVec.begin(), BinVec.end(), 0);
				});

			muios.lock();
			std::cout << "Thread " << ThreadNo + 1 << ". Finished reference unit " << ReferenceUnit + 1
				<< " vs target unit " << TargetUnit + 1<< ".\n";
			muios.unlock();
		}

	}

}

void Statistician::CloseFiles()
{
	if (CorrFile.bad())
		std::cout << "bad";

	else if (CorrFile.eof())
		std::cout << "eof";

	else if (CorrFile.fail())
		std::cout << "other fail";

	else if (CorrFile.good())
	{
		CorrFile.close();
		muios.lock();
		std::cout << "Output file was closed successfully\n";
		muios.unlock();
	}

	if (CorrFile.rdstate() == (std::ios_base::failbit | std::ios_base::eofbit))
	{
		muios.lock();
		std::cout << "stream state is eofbit\n";
		muios.unlock();
	}

}


//Func previous to Multithreading implementation. Rehabilitated only for debugging puroposes. Multithreading should be use for the actual analyses.
void Statistician::MasterSpikeCrossCorrSingleThread(int Stimulus, int ResampledSets, uint8_t ResamplingMethod, uint8_t StatTest, double PVal, bool ExcZeroLag)
{

	//NOTE: Im not convienced that STD matrices are the best way to deal with the problem. They are well allocated but anyway they may impact the performance,
	//STD problem may be solved with the use of other statistics instead of Z test (Fujisawa, 2008).
	//it is a posibility to implement fujisawa statistics but I need to try them first on MATLAB.

	//mu.lock();

	//Locked code to access common memory between threads
	int UnitsRef = OdorEx.GetUnitsRef();
	int UnitsTar = OdorEx.GetUnitsTar();
	int Trials = OdorEx.GetTrials();
	auto SLSRB = StimLockedSpikesRef.cbegin();
	auto SLSTB = StimLockedSpikesTar.cbegin();

	std::cout << "Stimulus: " << Stimulus + 1 << ", Ref: " << UnitsRef
		<< ", Tar: " << UnitsTar << ", Trials: " << Trials << "\n";

	//Put this thread to sleep just for debugging puposes.
	//std::this_thread::sleep_for(std::chrono::milliseconds(1000));
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<unsigned int> SpikesCountCorr(NoBins); //Main raw correlation Vector

	std::vector<std::vector<unsigned int>> SpikesCountResampled(ResampledSets, std::vector<unsigned int>(NoBins)); // Good! Resampling Matrix, this is annoying but necessary to obtain the standard deviation.
	std::vector<unsigned int> SpikesSTDCount(ResampledSets);

	//Vars when working with PermTest comp fujisawa, 2008.
	std::vector<uint32_t> LPWBand(NoBins); //
	std::vector<uint32_t> UPWBand(NoBins);
	int PValPlace = (int)std::ceil(double(ResampledSets * PVal) / 2.0);
	std::vector<std::vector<uint32_t>> LPWBands(PValPlace, std::vector<uint32_t>(NoBins));
	std::vector<std::vector<uint32_t>> UPWBands(PValPlace, std::vector<uint32_t>(NoBins));
	std::pair<uint32_t, uint32_t> GlobalBands(0, 0);
	//mu.unlock();

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//File for storing the sig data.
	std::ofstream CorrFile("Stimulus" + std::to_string(Stimulus + 1) + ".txt");
	std::ofstream JitteredMatrixFile("JitteredMatrix" + std::to_string(Stimulus + 1) + ".txt");

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
		if (ReferenceUnit == 145)
		{
			//Stimulus locked target spike train loop
			unsigned short TargetUnit = 1;
			for (auto TarTrain = SLSTB + ((__int64)Stimulus * UnitsTar),
				endTT = TarTrain + UnitsTar;
				TarTrain < endTT
				; ++TarTrain, TargetUnit++)
			{
				if (TargetUnit == 181)
				//if (ReferenceUnit < TargetUnit)
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
							switch (ResamplingMethod)
							{
							case SHUFFLING:
								SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, CountCorr); //Compute Corr.
								SpikeTrainShuffle(*RefTrialTrain, *TarTrialTrain, SpikesCountResampled, CountRes); //Compute Corr shuffling method.
								break;

							case BASICJITTER:
								SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, CountCorr);
								SpikeTrainBasicJitter(*RefTrialTrain, *TarTrialTrain, SpikesCountResampled, CountRes); //Compute Basic Jittering method.
								break;

							case INTERJITTER:
								SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, CountCorr);
								SpikeTrainIntervalJitter(SpikesCountCorr, SpikesCountResampled, CountRes);
								//SpikeTrainIntervalJitter4(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, SpikesCountResampled, CountCorr);
								//WriteToFileWorkerT(CorrFile, SpikesCountCorr); CorrFile << "\n";

								break;

							default:
								SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, CountCorr);
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

						for (auto BinVec = SpikesCountResampled.begin(), BinVecEnd = SpikesCountResampled.end();
							BinVec < BinVecEnd;
							++BinVec)
						{
							/*if (LPWVecit == LPWBands.cbegin())
							{
								WriteToFileWorkerT(JitteredMatrixFile, *BinVec);
								JitteredMatrixFile << "\n";
							}*/
							for (auto LPWit = LPWVecit->cbegin(), UPWit = UPWVecit->cbegin(), End = LPWVecit->cend(), ResDatait = BinVec->cbegin();
								LPWit < End; ++LPWit, ++UPWit, ++ResDatait)
							{
								if (*ResDatait < *LPWit || *ResDatait >* UPWit)
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
							break;
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


					/*WriteToFileWorkerT(JitteredMatrixFile, LPWBand);
					JitteredMatrixFile << "\n";
					WriteToFileWorkerT(JitteredMatrixFile, UPWBand);
					JitteredMatrixFile << "\n";
					JitteredMatrixFile << GlobalBands.first<< "\n";
					JitteredMatrixFile << GlobalBands.second << "\n";*/

					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					/////////////////////////////////////////////Permutation Test///////////////////////////////////////////////////////////////////////////////////

					if ((CountCorr != 0) && GoodData)
					{
						//Writing Significant Correlations to File.

						bool SigArray[4] = { false };

						GetSignifcantCorr(SpikesCountCorr, SigArray, GlobalBands, LPWBand, UPWBand);

						//mu.lock();
						CorrFile << 8 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
						WriteToFileWorkerT(CorrFile, SpikesCountCorr);
						WriteToFileWorkerT(CorrFile, LPWBand);
						WriteToFileWorkerT(CorrFile, UPWBand);
						CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";

						if ((SigArray[2] || SigArray[3]) && (SigArray[0] || SigArray[1]))
						{
							CorrFile << 7 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
							WriteToFileWorkerT(CorrFile, SpikesCountCorr);
							WriteToFileWorkerT(CorrFile, LPWBand);
							WriteToFileWorkerT(CorrFile, UPWBand);
							CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
						}
						else
						{
							if (SigArray[0] && SigArray[1])
							{
								//Code to store in txt files.
								CorrFile << 1 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
								WriteToFileWorkerT(CorrFile, SpikesCountCorr);
								WriteToFileWorkerT(CorrFile, LPWBand);
								WriteToFileWorkerT(CorrFile, UPWBand);
								CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
							}
							else if (SigArray[0])
							{
								CorrFile << 2 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
								WriteToFileWorkerT(CorrFile, SpikesCountCorr);
								WriteToFileWorkerT(CorrFile, LPWBand);
								WriteToFileWorkerT(CorrFile, UPWBand);
								CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
							}
							else if (SigArray[1])
							{
								CorrFile << 3 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
								WriteToFileWorkerT(CorrFile, SpikesCountCorr);
								WriteToFileWorkerT(CorrFile, LPWBand);
								WriteToFileWorkerT(CorrFile, UPWBand);
								CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
							}

							if (SigArray[2] && SigArray[3])
							{
								CorrFile << 4 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
								WriteToFileWorkerT(CorrFile, SpikesCountCorr);
								WriteToFileWorkerT(CorrFile, LPWBand);
								WriteToFileWorkerT(CorrFile, UPWBand);
								CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
							}
							else if (SigArray[2])
							{
								CorrFile << 5 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
								WriteToFileWorkerT(CorrFile, SpikesCountCorr);
								WriteToFileWorkerT(CorrFile, LPWBand);
								WriteToFileWorkerT(CorrFile, UPWBand);
								CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
							}
							else if (SigArray[3])
							{
								CorrFile << 6 << ", " << ReferenceUnit << ", " << TargetUnit << ", ";
								WriteToFileWorkerT(CorrFile, SpikesCountCorr);
								WriteToFileWorkerT(CorrFile, LPWBand);
								WriteToFileWorkerT(CorrFile, UPWBand);
								CorrFile << GlobalBands.first << ", " << GlobalBands.second << ", " << CountCorr << ", " << GoodAlpha << ", " << "\n";
							}
						}
						//mu.unlock();
					}

					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

					//Reseting Count Vectors and Matrix;
					std::fill(SpikesCountCorr.begin(), SpikesCountCorr.end(), 0);

					std::for_each(SpikesCountResampled.begin(), SpikesCountResampled.end(),
						[](std::vector<unsigned int>& BinVec)
						{
							std::fill(BinVec.begin(), BinVec.end(), 0);
						});

					//mu.lock();
					std::cout << "Stimulus " << Stimulus + 1 << ". Finished reference unit " << ReferenceUnit << " vs target unit " << TargetUnit << ".\n";
					//mu.unlock();
				}
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
		//mu.lock();
		std::cout << "Output file was closed successfully\n";
		//mu.unlock();
	}

	if (CorrFile.rdstate() == (std::ios_base::failbit | std::ios_base::eofbit))
	{
		//mu.lock();
		std::cout << "stream state is eofbit\n";
		//mu.unlock();
	}
} 
