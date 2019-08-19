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
		// aritmethic with custom made functions. STA stands for Spike Train Analysis.

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

void Statistician::SpikeTrainIntervalJitter(const std::vector<uint32_t>& reference, const std::vector<uint32_t>& target, std::vector<unsigned int>& Spikes, std::vector<std::vector<unsigned int>>& SpikesMatrix,unsigned int& Count)
{
	//To produce a partition of t into k values :

	//	- Generate k - 1 uniformly distributed values in the range[0, t].

	//	- Sort them, and add 0 at the beginning and t at the end.

	//	- Use the adjacent differences as the partition.


	//Somehow this is failing. Im getting just a normal distribution when I need a uniform dist.
	// The normal dist is evident even in MATLAB. when t < k , it tends to be a normal dist.
	//I need to implement the real test with pseudorandom spike samples not counts but it may be really slow.

	long long CurrentBinF;
	long long CurrentBinL;
	auto LBit = target.begin();
	auto UBit = target.begin();

	int HelperSize = NoBins / 2 + 1;

	std::vector<uint32_t> HelperLag(HelperSize);
	std::vector<uint32_t> HelperLead(HelperSize);

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
		auto LagIt = HelperLag.begin() + 1, LeadIt = HelperLead.begin() + 1;
		const auto LagEnd = HelperLag.end(), LeadEnd = HelperLead.end();

		/////This loops are written this way to avoid counting zero lag correlations. They are implemented using pointer
		// aritmethic with custom made functions. STA stands for Spike Train Analysis.

		uint32_t BinCount;

		for (; Bin < LastBin - (NoBins / 2) - 1; ++Bin, ++LagIt)
		{
			BinCount = (uint32_t)std::distance(First, Last);
			*LagIt = BinCount;
			*Bin += BinCount;

			CurrentBinF = CurrentBinL;
			CurrentBinL = CurrentBinF + BinSize;

			First = Last;
			STAUpperBoundTExc(Last, UBit, CurrentBinL);
		}

		BinCount = (uint32_t)std::distance(First, Last);
		*LagIt = BinCount;
		*Bin += BinCount;
		CurrentBinF = CurrentBinL;
		CurrentBinL = CurrentBinF + BinSize;

		STALowerBoundTExc(First, UBit, CurrentBinF);
		STAUpperBoundT(Last, UBit, CurrentBinL);

		++Bin;

		for (; Bin < LastBin; ++Bin, ++LeadIt)
		{
			BinCount = (uint32_t)std::distance(First, Last);
			*LeadIt = BinCount;
			*Bin += BinCount;
			CurrentBinF = CurrentBinL;
			CurrentBinL = CurrentBinF + BinSize;

			First = Last;
			STAUpperBoundT(Last, UBit, CurrentBinL);
		}

		LagIt = HelperLag.begin() + 1;
		LeadIt = HelperLead.begin() + 1;

		uint32_t LagCount = std::accumulate(LagIt, LagEnd,0);
		uint32_t LeadCount = std::accumulate(LeadIt, LeadEnd, 0);

		std::uniform_int_distribution<uint32_t> LagDist(0, LagCount);
		std::uniform_int_distribution<uint32_t> LeadDist(0, LeadCount);
		*(HelperLag.end() - 1) = LagCount;
		*(HelperLead.end() - 1) = LeadCount;


		for (auto JittVec = SpikesMatrix.begin(), SMEnd = SpikesMatrix.end(); JittVec < SMEnd; ++JittVec)
		{

			for (; LagIt < (LagEnd - 1); ++LagIt)
			{
				*LagIt = LagDist(Generator);
				//std::cout << *LagIt << "\n";
			}

			for (; LeadIt < (LeadEnd - 1); ++LeadIt)
			{
				*LeadIt = LeadDist(Generator);
			}

			std::sort(HelperLag.begin(), HelperLag.end());
			std::sort(HelperLead.begin(), HelperLead.end());

			LagIt = HelperLag.begin() + 1;
			LeadIt = HelperLead.begin() + 1;

			for (auto JittLagIt = JittVec->begin(), JittLagEnd = JittVec->begin() + (NoBins / 2); JittLagIt < JittLagEnd; ++JittLagIt, ++LagIt)
			{
				*JittLagIt += *LagIt - *(LagIt - 1);
			}

			for (auto JittLeadIt = JittVec->begin() + (NoBins / 2), JittLeadEnd = JittVec->end(); JittLeadIt < JittLeadEnd; ++JittLeadIt, ++LeadIt)
			{
				*JittLeadIt += *LeadIt - *(LeadIt - 1);
			}

			LagIt = HelperLag.begin() + 1;
			LeadIt = HelperLead.begin() + 1;
		}
	}
	Count += (unsigned int)reference.size();
}

void Statistician::SpikeTrainIntervalJitter2(std::vector<unsigned int>& Spikes, std::vector<std::vector<unsigned int>>& SpikesMatrix)
{

	int HelperSize = NoBins / 2 + 1;

	std::vector<uint32_t> HelperLag(HelperSize);
	std::vector<uint32_t> HelperLead(HelperSize);

	auto LagIt = HelperLag.begin() + 1;
	auto LeadIt = HelperLead.begin() + 1;
	auto LagEnd = HelperLag.end();
	auto LeadEnd = HelperLead.end();

	uint32_t LagCount = std::accumulate(Spikes.begin(), Spikes.begin() + NoBins / 2, 0);
	uint32_t LeadCount = std::accumulate(Spikes.begin() + NoBins / 2, Spikes.end(), 0);

	std::uniform_int_distribution<uint32_t> LagDist(0, LagCount);
	std::uniform_int_distribution<uint32_t> LeadDist(0, LeadCount);
	*(HelperLag.end() - 1) = LagCount;
	*(HelperLead.end() - 1) = LeadCount;

	for (auto JittVec = SpikesMatrix.begin(), SMEnd = SpikesMatrix.end(); JittVec < SMEnd; ++JittVec)
	{

		for (; LagIt < (LagEnd - 1); ++LagIt)
		{
			*LagIt = LagDist(Generator);
		}

		for (; LeadIt < (LeadEnd - 1); ++LeadIt)
		{
			*LeadIt = LeadDist(Generator);
		}

		std::sort(HelperLag.begin(), HelperLag.end());
		std::sort(HelperLead.begin(), HelperLead.end());

		LagIt = HelperLag.begin() + 1;
		LeadIt = HelperLead.begin() + 1;

		for (auto JittLagIt = JittVec->begin(), JittLagEnd = JittVec->begin() + (NoBins / 2); JittLagIt < JittLagEnd; ++JittLagIt, ++LagIt)
		{
			*JittLagIt += *LagIt - *(LagIt - 1);
		}

		for (auto JittLeadIt = JittVec->begin() + (NoBins / 2), JittLeadEnd = JittVec->end(); JittLeadIt < JittLeadEnd; ++JittLeadIt, ++LeadIt)
		{
			*JittLeadIt += *LeadIt - *(LeadIt - 1);
		}

		LagIt = HelperLag.begin() + 1;
		LeadIt = HelperLead.begin() + 1;
	}


}

void Statistician::SpikeTrainIntervalJitter3(std::vector<unsigned int>& Spikes, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count)
{
	uint32_t JitterInterval = 150; //This is gonna be 150 because I want an interval of 5 ms.
	uint32_t JitterCounts = (Epoch / JitterInterval) * 2;


	uint32_t LagCount = std::accumulate(Spikes.begin(), Spikes.begin() + NoBins / 2, 0);
	uint32_t LeadCount = std::accumulate(Spikes.begin() + NoBins / 2, Spikes.end(), 0);

	std::vector<uint32_t> JitterIntC(JitterCounts);

	//Fake RefSample is any element of the set {x:x>Epoch};
	uint32_t FakeRefSample = 100000;
	std::vector<uint32_t> FakeReference(1, FakeRefSample);
	std::vector<uint32_t> FakeTarget((long long)LagCount + (long long)LeadCount);

	std::uniform_int_distribution<uint32_t> LagDist(FakeRefSample - Epoch, FakeRefSample - 1);
	std::uniform_int_distribution<uint32_t> LeadDist(FakeRefSample + 1, FakeRefSample + Epoch);




	for (auto JittVec = SpikesMatrix.begin(), SMEnd = SpikesMatrix.end(); JittVec < SMEnd; ++JittVec)
	{
		auto FTIt = FakeTarget.begin();

		for (auto FTLagEnd = FakeTarget.begin() + LagCount; FTIt < FTLagEnd; ++FTIt)
		{
			*FTIt = LagDist(Generator);
			//std::cout << *FTIt << "\n";
		}

		for (auto FTLeadEnd = FakeTarget.end(); FTIt < FTLeadEnd; ++FTIt)
		{
			*FTIt = LeadDist(Generator);
		}
		std::sort(FakeTarget.begin(), FakeTarget.end());
		SpikeTrainCorr(FakeReference, FakeTarget, *JittVec, Count);
	}
}

void Statistician::SpikeTrainIntervalJitter4(const std::vector<uint32_t>& reference, const std::vector<uint32_t>& target, std::vector<unsigned int>& Spikes, std::vector<std::vector<unsigned int>>& SpikesMatrix, unsigned int& Count)
{

	long long CurrentBinF;
	long long CurrentBinL;
	auto LBit = target.begin();
	auto UBit = target.begin();

	int HelperSize = NoBins / 2;

	std::vector<uint32_t> HelperLag(HelperSize);
	std::vector<uint32_t> HelperLead(HelperSize);

	//std::cout << "Shuff: " << "\n";
	//for(auto Spike = reference.begin(), LastSpike = reference.end(); Spike < LastSpike; ++Spike)
	for (const uint32_t& Spike : reference)
	{
		CurrentBinF = (long long)Spike - Epoch; // Set the current bins for the lambda function.
		CurrentBinL = CurrentBinF + BinSize;

		std::uniform_int_distribution<uint32_t> LagDist(Spike - Epoch, Spike - 1);
		std::uniform_int_distribution<uint32_t> LeadDist(Spike + 1, (long long)Spike + (long long)Epoch);


		//Boundaries of the target spikes.

		STALowerBoundT(LBit, target.end(), CurrentBinF);
		STAUpperBoundT(UBit, target.end(), Spike + Epoch);

		auto First = LBit;
		auto Last = LBit;
		STAUpperBoundTExc(Last, UBit, CurrentBinL);

		//Ierators for Count Corr vec
		auto Bin = Spikes.begin(), LastBin = Spikes.end();
		auto LagIt = HelperLag.begin(), LeadIt = HelperLead.begin();
		const auto LagEnd = HelperLag.end(), LeadEnd = HelperLead.end();


		/////This loops are written this way to avoid counting zero lag correlations. They are implemented using pointer
		// aritmethic with custom made functions. STA stands for Spike Train Analysis.

		uint32_t BinCount;

		for (; Bin < LastBin - (NoBins / 2) - 1; ++Bin, ++LagIt)
		{
			BinCount = (uint32_t)std::distance(First, Last);
			*LagIt = BinCount;
			*Bin += BinCount;

			CurrentBinF = CurrentBinL;
			CurrentBinL = CurrentBinF + BinSize;

			First = Last;
			STAUpperBoundTExc(Last, UBit, CurrentBinL);
		}

		BinCount = (uint32_t)std::distance(First, Last);
		*LagIt = BinCount;
		*Bin += BinCount;
		CurrentBinF = CurrentBinL;
		CurrentBinL = CurrentBinF + BinSize;

		STALowerBoundTExc(First, UBit, CurrentBinF);
		STAUpperBoundT(Last, UBit, CurrentBinL);

		++Bin;

		for (; Bin < LastBin; ++Bin, ++LeadIt)
		{
			BinCount = (uint32_t)std::distance(First, Last);
			*LeadIt = BinCount;
			*Bin += BinCount;
			CurrentBinF = CurrentBinL;
			CurrentBinL = CurrentBinF + BinSize;

			First = Last;
			STAUpperBoundT(Last, UBit, CurrentBinL);
		}

		LagIt = HelperLag.begin(), LeadIt = HelperLead.begin();
		uint32_t LagCount = std::accumulate(LagIt, LagEnd, 0);
		uint32_t LeadCount = std::accumulate(LeadIt, LeadEnd, 0);

		//Fake RefSample is any element of the set {x:x>Epoch};
		std::vector<uint32_t> FakeReference(1, Spike);
		std::vector<uint32_t> FakeTarget((long long)LagCount + (long long)LeadCount);

		for (auto JittVec = SpikesMatrix.begin(), SMEnd = SpikesMatrix.end(); JittVec < SMEnd; ++JittVec)
		{
			auto FTIt = FakeTarget.begin();


			std::vector<uint32_t> Checker((long long)LagCount + (long long)LeadCount + 1L);
			auto ChIt = Checker.begin();

			for (auto FTLagEnd = FakeTarget.begin() + LagCount; FTIt < FTLagEnd; ++FTIt)
			{
				uint32_t LagSample = LagDist(Generator);

				while (std::any_of(Checker.cbegin(), Checker.cend(), [&LagSample](const uint32_t& Ch) {return Ch == LagSample; }))
				{
					LagSample = LagDist(Generator);
				}

				*FTIt = LagSample;
				*ChIt = LagSample;
				ChIt++;
			}

			for (auto FTLeadEnd = FakeTarget.end(); FTIt < FTLeadEnd; ++FTIt)
			{
				//*FTIt = LeadDist(Generator);
				uint32_t LeadSample = LeadDist(Generator);

				while (std::any_of(Checker.cbegin(), Checker.cend(), [&LeadSample](const uint32_t& Ch) {return Ch == LeadSample; }))
				{
					LeadSample = LeadDist(Generator);
				}

				*FTIt = LeadSample;
				*ChIt = LeadSample;
				ChIt++;
			}
			std::sort(FakeTarget.begin(), FakeTarget.end());
			unsigned int FakeCount = 0;
			SpikeTrainCorr(FakeReference, FakeTarget, *JittVec, FakeCount);
		}
	}
	Count += (unsigned int)reference.size();
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
	int PValPlace = (int)std::ceil(double(ResampledSets * PVal) / 2.0);
	std::vector<std::vector<uint32_t>> LPWBands(PValPlace, std::vector<uint32_t>(NoBins));
	std::vector<std::vector<uint32_t>> UPWBands(PValPlace, std::vector<uint32_t>(NoBins));
	std::pair<uint32_t, uint32_t> GlobalBands(0, 0);
	mu.unlock();

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
		if (ReferenceUnit == 37)
		{
			//Stimulus locked target spike train loop
			unsigned short TargetUnit = 1;
			for (auto TarTrain = SLSTB + ((__int64)Stimulus * UnitsTar),
				endTT = TarTrain + UnitsTar;
				TarTrain < endTT
				; ++TarTrain, TargetUnit++)
			{
				if (TargetUnit == 49)
				//if(ReferenceUnit != TargetUnit)
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
								//SpikeTrainCorr(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, CountCorr);
								//SpikeTrainIntervalJitter3(SpikesCountCorr, SpikesCountResampled,CountRes);
								SpikeTrainIntervalJitter4(*RefTrialTrain, *TarTrialTrain, SpikesCountCorr, SpikesCountResampled, CountCorr);
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
							if (LPWVecit == LPWBands.cbegin())
							{ 
								//WriteToFileWorkerT(JitteredMatrixFile, *BinVec);
								//JitteredMatrixFile << "\n";
							}
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

						mu.lock();
						
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
		std::cout << "Output file was closed successfully\n";
		mu.unlock();
	}

	if (CorrFile.rdstate() == (std::ios_base::failbit | std::ios_base::eofbit))
	{
		mu.lock();
		std::cout << "stream state is eofbit\n";
		mu.unlock();
	}
}