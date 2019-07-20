#include "Statistician.h"
#include <algorithm>
#include <future>
#include <iostream>
#include <numeric>
#include <cmath>


constexpr char SHUFFLING = 0;
constexpr char JITTERING = 1;

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

void Statistician::SpikeTrainJitter(const std::vector<double>& reference, std::vector<double> target, std::vector<std::vector<unsigned int>>& Spikes, unsigned int& Count)
{
}

void Statistician::SpikeTrainShuffle(const std::vector<double>& reference, std::vector<double> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count)
{

	double TargetMin;
	double TargetMax;
	double TargetRandom;

	//Setting Pseudo Random Number Uniform Distribution
	std::uniform_int_distribution<int> distribution(1, target.size() - 1);

	for (auto Spikes = SpikesMatrix.begin(), SMEnd = SpikesMatrix.end(); Spikes < SMEnd ; ++Spikes)
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
			*Spikes,
			Count);
	}
}

void Statistician::MasterSpikeCrossCorr()
{

	//Nested loops for running the whole analysis. There may be some improvement specially in the las loop if the data is parsed better from matlab.
	//first loop might be implemented in a multithreading design.

	for (int Stimulus = 0; Stimulus < OdorEx.GetStimuli(); Stimulus++)
	{
		for (auto RefTrain = StimLockedSpikesRef.cbegin() + Stimulus * OdorEx.GetUnitsRef(),
			endRT = RefTrain + OdorEx.GetUnitsRef();
			RefTrain < endRT
			; ++RefTrain)
		{
			for (auto TarTrain = StimLockedSpikesTar.cbegin() + Stimulus * OdorEx.GetUnitsTar(),
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
		*CurrentThread = std::async(std::launch::async, &Statistician::MasterSpikeCrossCorrWorker, this, Stimulus, );
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

void Statistician::MasterSpikeCrossCorrWorker(int Stimulus, int ResampledSets, char ResamplingMethod)
{

	//NOTE: Im not convienced that STD matrices are the best way to deal with the problem. They are well allocated but anyway they may impact the performance,
	//STD problem may be solved with the use of other statistics instead of Z test (Fujisawa, 2008).
	//it is a posibility to implement fujisawa statistics but I need to try them first on MATLAB.


	mu.lock();
	//Locked code to access common memory between threads.
	std::cout << "adq vars valve :" << Stimulus << std::endl;

	int UnitsRef = OdorEx.GetUnitsRef();
	int UnitsTar = OdorEx.GetUnitsTar();
	int Trials = OdorEx.GetTrials();
	auto SLSRB = StimLockedSpikesRef.cbegin();
	auto SLSTB = StimLockedSpikesTar.cbegin();

	std::cout << "end adq vars valve :" << Stimulus << std::endl;

	//Put this thread to sleep just for debugging puposes.
	std::this_thread::sleep_for(std::chrono::milliseconds(1000));

	mu.unlock();

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<unsigned int> SpikesCountCorr(NoBins); //Main raw correlation Vector
	std::vector<double> SpikesPCorr(NoBins); // Vector for probabilities and Z scores. P stands for probability.

	std::vector<std::vector<unsigned int>> SpikesCountResampled(ResampledSets,std::vector<unsigned int>(NoBins)); // Good! Resampling Matrix, this is annoying but necessary to obtain the standard deviation.
	std::vector<unsigned int> SpikesSTDCount(ResampledSets);
	std::vector<double> SpikesSTDResampled(NoBins); // STD vector. STD is obtained across ResampledSets of mean trials.
	std::vector<double> SpikesPResampled(NoBins); // Vector for probabilities scores.

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	//Nested loops for running the whole analysis. There may be some improvement specially in the las loop if the data is parsed better from matlab.

	//Stimulus locked reference spike train loop
	for (auto RefTrain = SLSRB + Stimulus * UnitsRef,
		endRT = RefTrain + UnitsRef;
		RefTrain < endRT
		; ++RefTrain)
	{
		//Stimulus locked target spike train loop
		for (auto TarTrain = SLSTB + Stimulus * UnitsTar,
			endTT = TarTrain + UnitsTar;
			TarTrain < endTT
			; ++TarTrain)
		{
			auto RefTrialTrain = RefTrain; //this is the downside of the way I parse the matlab data.
			auto TarTrialTrain = TarTrain; //Aux vars to prevent modification of original vars.
			unsigned int CountCorr = 0;
			unsigned int CountRes = 0;
			bool GoodResampling = true;

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


			
			//Mean and STD of the matrix and vectors of the choosen resampling method.///////////
			CountRes /= ResampledSets;
			for (int Bin = 0; Bin < NoBins; Bin++)
			{
				auto STDCount = SpikesSTDCount.begin();
				auto STDCountEnd = SpikesSTDCount.end();

				for (auto BinVec = SpikesCountResampled.cbegin(), BinVecEnd = SpikesCountResampled.cend();
					BinVec < BinVecEnd;
					++BinVec, ++STDCount)
				{
					*STDCount = *(BinVec->begin() + Bin);
				}

				//Math for the params that are needed by the Z Test
				double BinMean = (double)std::accumulate(STDCount, STDCountEnd, 0) / (double)ResampledSets;

				//Necesary check for unpopulated resampled correlograms. false positives can be assumed if this is not checked, although this is not the best way to code it. Bad design.
				if (BinMean == 0)
					GoodResampling = false; break;

				double BinVariance = 0.0;

				for (STDCount = SpikesSTDCount.begin(); STDCount < STDCountEnd; ++STDCount)
				{
					BinVariance += ((double)(*STDCount) - BinMean) * ((double)(*STDCount) - BinMean);
				}
				BinVariance /= (double)ResampledSets; // this is Variance over N. Matlab uses Bessels correction to compute STD.

				*(SpikesSTDResampled.begin() + Bin) = std::sqrt(BinVariance) / (double)CountRes; // Stand deviation to my STD.
				*(SpikesPResampled.begin() + Bin) = BinMean / (double)CountRes;

				STDCount = SpikesSTDCount.begin(); // Reseting the iterator of the vector.
			}
			/////////////////////////////////////////////////////////////////////////////////////


			if (GoodResampling && CountCorr != 0)
			{
				
				double MeanSTD = (double)std::accumulate(SpikesSTDResampled.begin(), SpikesSTDResampled.end(), 0.0) / (double)SpikesSTDResampled.size();
				
				//Fill the Probability Vector.
				std::transform(SpikesCountCorr.begin(), SpikesCountCorr.end(),
					SpikesPCorr.begin(),
					[CountCorr](unsigned int Bin) -> double { return (double)Bin / (double)CountCorr; });


				std::transform(SpikesPCorr.begin(), SpikesPCorr.begin(), SpikesPResampled.begin(), SpikesPCorr.begin(),
					[MeanSTD](double& PBin, double& MeanBin) {return (PBin - MeanBin) / MeanSTD; });




				//CHECAR TODOS LOS FOR EACH Y TRANSFORMS A LA DE YA!!!!!


			}
		}
	}
}