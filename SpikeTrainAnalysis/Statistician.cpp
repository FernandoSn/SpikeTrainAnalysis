#include "Statistician.h"
#include <algorithm>
#include <future>
#include <iostream>
#include <numeric>
#include <cmath>
#include <string>




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

Statistician::Statistician(std::string FileName, double Interval, int BinSize, int Epoch)
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

		std::cout << "Shuff: " << "\n";
	for(auto Spike = reference.begin(), LastSpike = reference.end(); Spike < LastSpike; ++Spike)
	//for (const double& Spike : reference)
	{
		CurrentBinF = *Spike - EpochSec; // Set the current bins for the lambda function.
		CurrentBinL = CurrentBinF + BinSizeSec;

		for (auto Bin = Spikes.begin(), LastBin = Spikes.end(); Bin < LastBin; ++Bin)
		//for (unsigned int& Bin : Spikes)
		{
			*Bin += (unsigned int)std::count_if(target.begin(), //std library stuff, very convenient and fast.
				target.end(),
				[&CurrentBinF, &CurrentBinL](double TargetSpike) 
				{
					return TargetSpike > CurrentBinF && TargetSpike <= CurrentBinL; 
				}
			);

			CurrentBinF = CurrentBinL;
			CurrentBinL = CurrentBinF + BinSizeSec;
		}
	}
		std::cout << "Shuff2: " << "\n";
	Count += (unsigned int)reference.size();
}

void Statistician::SpikeTrainJitter(const std::vector<double>& reference, std::vector<double> target, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count)
{

	//Setting Pseudo Random Number Uniform Distribution for jittering [-5ms, 5ms] (Fujisawa, 2018)
	std::uniform_real_distribution<double> distribution(-0.005, 0.005);

	for (auto Spikes = SpikesMatrix.begin(), SMEnd = SpikesMatrix.end(); Spikes < SMEnd; ++Spikes)
	{

		double Jitter = distribution(Generator);

		std::for_each(target.begin(),
			target.end(),
			[&Jitter](double& Spike) { Spike += Jitter; });

		SpikeTrainCorr(reference,
			target,
			*Spikes,
			Count);
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

void Statistician::MasterSpikeCrossCorrDeprecated(int ResampledSets, unsigned char ResamplingMethod, double ZThresh, bool ExcZeroLag)
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<unsigned int> SpikesCountCorr(NoBins); //Main raw correlation Vector
	std::vector<double> SpikesPCorr(NoBins); // Vector for probabilities and Z scores. P stands for probability.

	std::vector<std::vector<unsigned int>> SpikesCountResampled(ResampledSets, std::vector<unsigned int>(NoBins)); // Good! Resampling Matrix, this is annoying but necessary to obtain the standard deviation.
	std::vector<unsigned int> SpikesSTDCount(ResampledSets);
	std::vector<double> SpikesSTDResampled(NoBins); // STD vector. STD is obtained across ResampledSets of mean trials.
	std::vector<double> SpikesPResampled(NoBins); // Vector for probabilities scores.

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::ofstream CorrFile("TotalStimulus.SigCorr", std::ios::binary);

	//Check if we want to exclude "zero lag" correlations.
	int BinExcluded = 0;
	if (ExcZeroLag)
		BinExcluded = 1;

	//Nested loops for running the whole analysis. There may be some improvement specially in the las loop if the data is parsed better from matlab.
	//first loop might be implemented in a multithreading design.

	for (unsigned short Stimulus = 0; Stimulus < (unsigned short)OdorEx.GetStimuli(); Stimulus++)
	{
		unsigned short ReferenceUnit = 1;
		for (auto RefTrain = StimLockedSpikesRef.cbegin() + (long long)Stimulus * OdorEx.GetUnitsRef(),
			endRT = RefTrain + OdorEx.GetUnitsRef();
			RefTrain < endRT
			; ++RefTrain, ReferenceUnit++)
		{
			unsigned short TargetUnit = 1;
			for (auto TarTrain = StimLockedSpikesTar.cbegin() + (long long)Stimulus * OdorEx.GetUnitsRef(),
				endTT = TarTrain + OdorEx.GetUnitsTar();
				TarTrain < endTT
				; ++TarTrain, TargetUnit++)
			{
				auto RefTrialTrain = RefTrain; //this is the downside of the way I parse the matlab data.
				auto TarTrialTrain = TarTrain; //Aux vars to prevent modification of original vars.
				unsigned int CountCorr = 0;
				unsigned int CountRes = 0;
				bool GoodResampling = true;

				for (int Trial = 0; Trial < OdorEx.GetTrials(); 
					RefTrialTrain += OdorEx.GetUnitsRef(), TarTrialTrain += OdorEx.GetUnitsTar(), Trial++)
				{
					if ((RefTrialTrain->size() != 0 && TarTrialTrain->size() != 0))
					{
						SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, CountCorr); //Compute Corr.
						//std::cout << ".\n";
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

				//Mean and STD of the matrix and vectors of the choosen resampling method.///////////
				CountRes /= ResampledSets; //This needs to be divided into ResampledSets because that is the size of the Matrix, is not a vector anymore.
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
					BinVariance /= ((double)ResampledSets - 1.0); // this is Variance over N. Matlab uses Bessels correction to compute STD.

					*(SpikesSTDResampled.begin() + Bin) = std::sqrt(BinVariance) / (double)CountRes; // Stand deviation to my STD vector.
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
						[&CountCorr](unsigned int& Bin) -> double { return (double)Bin / (double)CountCorr; });

					//Z Transform.
					std::transform(SpikesPCorr.begin(), SpikesPCorr.end(), SpikesPResampled.begin(), SpikesPCorr.begin(),
						[&MeanSTD](double& PBin, double& MeanBin) -> double { return (PBin - MeanBin) / MeanSTD; });


					//Writing Significant Correlations to File.
					bool LeadEx = std::any_of(SpikesPCorr.end() - (SpikesPCorr.size() / 2) + BinExcluded, SpikesPCorr.end(),
						[&ZThresh](double& ZValue) {return ZValue > ZThresh; });
					bool LagEx = std::any_of(SpikesPCorr.begin(), SpikesPCorr.begin() + (SpikesPCorr.size() / 2) - BinExcluded,
						[&ZThresh](double& ZValue) {return ZValue > ZThresh; });
					bool LeadIn = std::any_of(SpikesPCorr.end() - (SpikesPCorr.size() / 2) + BinExcluded, SpikesPCorr.end(),
						[&ZThresh](double& ZValue) {return ZValue < -ZThresh; });
					bool LagIn = std::any_of(SpikesPCorr.begin(), SpikesPCorr.begin() + (SpikesPCorr.size() / 2) - BinExcluded,
						[&ZThresh](double& ZValue) {return ZValue < -ZThresh; });

					if (LeadEx && LagEx)
					{
						unsigned short CorrType = 1;
						CorrFile.write(reinterpret_cast<char*>(&Stimulus), 2);
						CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
						CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
						CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
					}
					else if (LeadEx)
					{
						unsigned short CorrType = 2;
						CorrFile.write(reinterpret_cast<char*>(&Stimulus), 2);
						CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
						CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
						CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
					}
					else if (LagEx)
					{
						unsigned short CorrType = 3;
						CorrFile.write(reinterpret_cast<char*>(&Stimulus), 2);
						CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
						CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
						CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
					}

					if (LeadIn && LagIn)
					{
						unsigned short CorrType = 4;
						CorrFile.write(reinterpret_cast<char*>(&Stimulus), 2);
						CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
						CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
						CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
					}
					else if (LeadIn)
					{
						unsigned short CorrType = 5;
						CorrFile.write(reinterpret_cast<char*>(&Stimulus), 2);
						CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
						CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
						CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
					}
					else if (LagIn)
					{
						unsigned short CorrType = 6;
						CorrFile.write(reinterpret_cast<char*>(&Stimulus), 2);
						CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
						CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
						CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
					}

					if ((LeadIn || LagIn) && (LeadEx || LagEx))
					{
						unsigned short CorrType = 7;
						CorrFile.write(reinterpret_cast<char*>(&Stimulus), 2);
						CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
						CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
						CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
					}
				}

				//Reseting Count Vectors and Matrix;
				std::fill(SpikesCountCorr.begin(), SpikesCountCorr.end(), 0);

				std::for_each(SpikesCountResampled.begin(), SpikesCountResampled.end(),
					[](std::vector<unsigned int>& BinVec)
					{
						std::fill(BinVec.begin(), BinVec.end(), 0);
					});

				std::cout << "Stimulus " << Stimulus << ". Finished reference unit " << ReferenceUnit << " vs target unit " << TargetUnit << ".\n";
			}
		}
	}
}

void Statistician::RunSingleThread(int ResampledSets, unsigned char ResamplingMethod, double ZThresh, bool ExcZeroLag)
{
	for (int Stimulus = 0; Stimulus < OdorEx.GetStimuli(); Stimulus++)
	{
		MasterSpikeCrossCorrWorker(Stimulus, ResampledSets, ResamplingMethod, ZThresh, ExcZeroLag);
	}
}

void Statistician::RunThreadPool(int ResampledSets, unsigned char ResamplingMethod, double ZThresh, bool ExcZeroLag)
{
	//ThreadPool with n threads. n = number of odors.
	std::vector<std::future<void>> ThreadPool(OdorEx.GetStimuli());

	auto CurrentThread = ThreadPool.begin();
	auto EndThread = ThreadPool.end();
	
	for (int Stimulus = 0; CurrentThread < EndThread; ++CurrentThread, Stimulus++)
	{
		*CurrentThread = std::async(std::launch::async, &Statistician::MasterSpikeCrossCorrWorker,
			this, Stimulus, ResampledSets, ResamplingMethod, ZThresh, ExcZeroLag);
	}

	while (true)
	{
		//This while loop freezes at the end, I need to implement in another way, but my thread pool works as the single thread.
		for (CurrentThread = ThreadPool.begin(); CurrentThread < EndThread; ++CurrentThread)
		{
			if (CurrentThread->valid())
				CurrentThread->get();
		}
	}
}

void Statistician::MasterSpikeCrossCorrWorker(int Stimulus, int ResampledSets, unsigned char ResamplingMethod, double ZThresh, bool ExcZeroLag)
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

	std::cout << "Stimulus: " << Stimulus << ", Ref: " << UnitsRef 
		<< ", Tar: " << UnitsTar << ", Trials: " << Trials << "\n";

	//Put this thread to sleep just for debugging puposes.
	//std::this_thread::sleep_for(std::chrono::milliseconds(1000));

	mu.unlock();

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::vector<unsigned int> SpikesCountCorr(NoBins); //Main raw correlation Vector
	std::vector<double> SpikesPCorr(NoBins); // Vector for probabilities and Z scores. P stands for probability.

	std::vector<std::vector<unsigned int>> SpikesCountResampled(ResampledSets,std::vector<unsigned int>(NoBins)); // Good! Resampling Matrix, this is annoying but necessary to obtain the standard deviation.
	std::vector<unsigned int> SpikesSTDCount(ResampledSets);
	std::vector<double> SpikesSTDResampled(NoBins); // STD vector. STD is obtained across ResampledSets of mean trials.
	std::vector<double> SpikesPResampled(NoBins); // Vector for probabilities scores.

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::ofstream CorrFile("Stimulus" + std::to_string(Stimulus + 1) + ".SigCorr", std::ios::binary);

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
			bool GoodResampling = true;

			mu.lock();
			std::cout << "Stimulus " << Stimulus + 1 << ". Finished reference unit " << ReferenceUnit << " vs target unit " << TargetUnit << ".\n";
			mu.unlock();

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
			CountRes /= ResampledSets; //This needs to be divided into ResampledSets because that is the size of the Matrix, is not a vector anymore.
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

				//Math for the params that are needed by the Z Test
				double BinMean = (double)std::accumulate(SpikesSTDCount.begin(), SpikesSTDCount.end(), 0) / (double)ResampledSets;

				//Necesary check for unpopulated resampled correlograms. false positives can be assumed if this is not checked, although this is not the best way to code it. Bad design.
				if (BinMean == 0)
					GoodResampling = false;

				double BinVariance = 0.0;

				for (STDCount = SpikesSTDCount.begin(); STDCount < STDCountEnd; ++STDCount)
				{
					BinVariance += ((double)(*STDCount) - BinMean) * ((double)(*STDCount) - BinMean);
				}
				BinVariance /= ((double)ResampledSets-1.0); // this is Variance over N. Matlab uses Bessels correction to compute STD.

				*(SpikesSTDResampled.begin() + Bin) = std::sqrt(BinVariance) / (double)CountRes; // Stand deviation to my STD vector.
				*(SpikesPResampled.begin() + Bin) = BinMean / (double)CountRes;

				STDCount = SpikesSTDCount.begin(); // Reseting the iterator of the vector.
			}
			/////////////////////////////////////////////////////////////////////////////////////

			if (GoodResampling && (CountCorr != 0))
			{
				
				double MeanSTD = (double)std::accumulate(SpikesSTDResampled.begin(), SpikesSTDResampled.end(), 0.0) / (double)SpikesSTDResampled.size();
				
				//Fill the Probability Vector.
				std::transform(SpikesCountCorr.begin(), SpikesCountCorr.end(),
					SpikesPCorr.begin(),
					[&CountCorr](unsigned int& Bin) -> double { return (double)Bin / (double)CountCorr; });

				//Z Transform.
				std::transform(SpikesPCorr.begin(), SpikesPCorr.end(), SpikesPResampled.begin(), SpikesPCorr.begin(),
					[&MeanSTD](double& PBin, double& MeanBin) -> double { return (PBin - MeanBin) / MeanSTD; });


				//Writing Significant Correlations to File.
				bool LeadEx = std::any_of(SpikesPCorr.end() - (SpikesPCorr.size() / 2) + BinExcluded, SpikesPCorr.end(),
					[&ZThresh](double& ZValue) {return ZValue > ZThresh; });
				bool LagEx = std::any_of(SpikesPCorr.begin(), SpikesPCorr.begin() + (SpikesPCorr.size() / 2) - BinExcluded,
					[&ZThresh](double& ZValue) {return ZValue > ZThresh; });
				bool LeadIn = std::any_of(SpikesPCorr.end() - (SpikesPCorr.size() / 2) + BinExcluded, SpikesPCorr.end(),
					[&ZThresh](double& ZValue) {return ZValue < -ZThresh; });
				bool LagIn = std::any_of(SpikesPCorr.begin(), SpikesPCorr.begin() + (SpikesPCorr.size() / 2) - BinExcluded,
					[&ZThresh](double& ZValue) {return ZValue < -ZThresh; });


				if (LeadEx && LagEx)
				{
					unsigned short CorrType = 1;
					CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
					CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
					CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
				}
				else if (LeadEx)
				{
					unsigned short CorrType = 2;
					CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
					CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
					CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
				}
				else if (LagEx)
				{
					unsigned short CorrType = 3;
					CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
					CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
					CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
				}

				if (LeadIn && LagIn)
				{
					unsigned short CorrType = 4;
					CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
					CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
					CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
				}
				else if (LeadIn)
				{
					unsigned short CorrType = 5;
					CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
					CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
					CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
				}
				else if (LagIn)
				{
					unsigned short CorrType = 6;
					CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
					CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
					CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
				}

				if ((LeadIn || LagIn) && (LeadEx || LagEx))
				{
					unsigned short CorrType = 7;
					CorrFile.write(reinterpret_cast<char*>(&CorrType), 2);
					CorrFile.write(reinterpret_cast<char*>(&ReferenceUnit), 2);
					CorrFile.write(reinterpret_cast<char*>(&TargetUnit), 2);
				}
			}

			//Reseting Count Vectors and Matrix;
			std::fill(SpikesCountCorr.begin(), SpikesCountCorr.end(), 0);

			std::for_each(SpikesCountResampled.begin(), SpikesCountResampled.end(),
				[](std::vector<unsigned int>& BinVec)
				{
					std::fill(BinVec.begin(), BinVec.end(), 0);
				});

			/*mu.lock();
			std::cout << "Stimulus " << Stimulus + 1 <<". Finished reference unit "<< ReferenceUnit << " vs target unit " << TargetUnit <<  ".\n";
			mu.unlock();*/
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
		std::cout << "Data file was closed successfully";
		std::cin.get();
		mu.unlock();
	}

	if (CorrFile.rdstate() == (std::ios_base::failbit | std::ios_base::eofbit))
	{
		mu.lock();
		std::cout << "stream state is eofbit\n";
		mu.unlock();
	}


}
