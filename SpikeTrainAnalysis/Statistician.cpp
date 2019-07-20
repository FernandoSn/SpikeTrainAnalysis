#include "Statistician.h"
#include <algorithm>
#include <future>
#include <iostream>


Statistician::Statistician(std::string FileName, int BinSize, int Epoch)
	:
	BinSize(BinSize),
	Epoch(Epoch),
	NoBins((Epoch / BinSize) * 2),
	BinSizeSec(BinSize / 1000.0),
	EpochSec(Epoch / 1000.0),
	OdorEx(FileName),
	Reference(OdorEx.RDataFile(), OdorEx.GetUnitsRef(), OdorEx.GetRefSizePos(), OdorEx.GetRefTrainPos()),
	Target(OdorEx.RDataFile(), OdorEx.GetUnitsTar(), OdorEx.GetTarSizePos(), OdorEx.GetTarTrainPos()),
	StimLockedSpikesRef((long long)OdorEx.GetStimuli() * (long long)OdorEx.GetMagnitudes() * (long long)OdorEx.GetTrials() * (long long)OdorEx.GetUnitsRef()),
	StimLockedSpikesTar((long long)OdorEx.GetStimuli() * (long long)OdorEx.GetMagnitudes() * (long long)OdorEx.GetTrials() * (long long)OdorEx.GetUnitsTar()),
	Generator(Rd())
{
	SetStimLockedSpikes();
}

Statistician::Statistician(std::string FileName, double Interval, int BinSize, int Epoch)
	:
	BinSize(BinSize),
	Epoch(Epoch),
	NoBins((Epoch / BinSize) * 2),
	BinSizeSec(BinSize / 1000.0),
	EpochSec(Epoch / 1000.0),
	OdorEx(FileName),
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

	for (const double& Spike : reference)
	{
		CurrentBinF = Spike - EpochSec; // Set the current bins for the lambda function.
		CurrentBinL = CurrentBinF + BinSizeSec;

		for (unsigned int& Bin : Spikes)
		{
			Bin += std::count_if(target.begin(), //std library stuff, very convenient and fast.
				target.end(),
				[&CurrentBinF, &CurrentBinL](double TargetSpike) 
				{
					return TargetSpike > CurrentBinF && TargetSpike <= CurrentBinL; 
				}
			);

			CurrentBinF = CurrentBinL;
			CurrentBinL = +BinSizeSec;
		}
	}
	Count += reference.size();
}

void Statistician::SpikeTrainShuffle(const std::vector<double>& reference, std::vector<double> target, std::vector<unsigned int>& Spikes, unsigned int& Count)
{

	double TargetMin;
	double TargetMax;
	double TargetRandom;

	//Setting Pseudo Random Number Uniform Distribution
	std::uniform_int_distribution<int> distribution(1, target.size() - 1);

	for (int i = 0; i < Shuffles; i++)
	{
		
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
			Spikes,
			Count);

	}
}

void Statistician::MasterSpikeCrossCorr()
{

	//Nested loops for running the whole analysis. There may be some improvement specially in the las loop if the data is parsed better from matlab.
	//first loop might be implemented in a multithreading design.

	for (int Stimulus = 0; Stimulus < OdorEx.GetStimuli(); Stimulus++)
	{
		for (auto RefTrain = StimLockedSpikesRef.begin() + Stimulus * OdorEx.GetUnitsRef(),
			endRT = RefTrain + OdorEx.GetUnitsRef();
			RefTrain < endRT
			; ++RefTrain)
		{
			for (auto TarTrain = StimLockedSpikesTar.begin() + Stimulus * OdorEx.GetUnitsTar(),
				endTT = TarTrain + OdorEx.GetUnitsTar();
				TarTrain < endTT
				; ++TarTrain)
			{
				auto RefTrialTrain = RefTrain;
				auto TarTrialTrain = TarTrain;

				for (int Trial = 0; Trial < OdorEx.GetTrials(); 
					RefTrialTrain += OdorEx.GetUnitsRef(), TarTrialTrain += OdorEx.GetUnitsTar(), Trial++)
				{
					if ((RefTrialTrain->size() != 0 && TarTrialTrain->size() != 0))
					{










					}
				}
			}
		}
	}
}

void Statistician::InitInterns()
{
	//ThreadPool with n threads. n = number of odors.
	std::vector<std::future<void>> ThreadPool(OdorEx.GetStimuli());

	auto CurrentThread = ThreadPool.begin();
	auto EndThread = ThreadPool.end();
	
	for (int Stimulus = 0; CurrentThread < EndThread; ++CurrentThread, Stimulus++)
	{
		*CurrentThread = std::async(std::launch::async, &Statistician::MasterSpikeCrossCorrWorker, this, Stimulus);
	}

	while (true)
	{
		for (CurrentThread = ThreadPool.begin(); CurrentThread < EndThread; ++CurrentThread)
		{
			if (CurrentThread->valid())
				CurrentThread->get();
		}
	}
}

void Statistician::MasterSpikeCrossCorrWorker(int Stimulus)
{

	mu.lock();

	std::cout << "adq vars valve :" << Stimulus << std::endl;

	int UnitsRef = OdorEx.GetUnitsRef();
	int UnitsTar = OdorEx.GetUnitsTar();
	int Trials = OdorEx.GetTrials();
	auto SLSRB = StimLockedSpikesRef.begin();
	auto SLSTB = StimLockedSpikesTar.begin();

	std::cout << "end adq vars valve :" << Stimulus << std::endl;

	std::this_thread::sleep_for(std::chrono::milliseconds(1000));


	mu.unlock();


	std::vector<unsigned int> SpikesCountCorr(NoBins);
	std::vector<unsigned int> SpikesCountShift(NoBins);
	std::vector<unsigned int> SpikesCountShuffle(NoBins);
	std::vector<unsigned int> SpikesCountJitter(NoBins);

	std::vector<double> SpikesPCorr(NoBins);
	std::vector<double> SpikesPShift(NoBins);
	std::vector<double> SpikesPShuffle(NoBins);
	std::vector<double> SpikesPJitter(NoBins);

	struct
	{
		unsigned int Corr = 0;
		unsigned int Jitter = 0;
		unsigned int Shift = 0;
		unsigned int Shuffle = 0;
	}Counts;

	//Nested loops for running the whole analysis. There may be some improvement specially in the las loop if the data is parsed better from matlab.
	for (auto RefTrain = SLSRB + Stimulus * UnitsRef,
		endRT = RefTrain + UnitsRef;
		RefTrain < endRT
		; ++RefTrain)
	{
		for (auto TarTrain = SLSTB + Stimulus * UnitsTar,
			endTT = TarTrain + UnitsTar;
			TarTrain < endTT
			; ++TarTrain)
		{
			auto RefTrialTrain = RefTrain;
			auto TarTrialTrain = TarTrain;

			for (int Trial = 0; Trial < Trials;
				RefTrialTrain += UnitsRef, TarTrialTrain += UnitsTar, Trial++)
			{
				if ((RefTrialTrain->size() != 0 && TarTrialTrain->size() != 0))
				{

					SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, Counts.Corr);

					SpikeTrainShuffle(*RefTrialTrain, *TarTrialTrain, SpikesCountShuffle, Counts.Shuffle);

				}
			}








		}
	}
}