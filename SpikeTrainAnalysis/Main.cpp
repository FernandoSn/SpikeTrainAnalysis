#include <iostream>
#include <chrono>
#include <thread>
#include <string>
#include <fstream>

#include "Statistician.h"
#include <numeric>
#include <iterator>
int main()
{
	//This is not professional code. to my taste, the classes are well designed maybe the thread pool or MasterCorr on Statistician
	//needs more work because a have a ton of code in one method but its acceptable.
	//The part that is really bad is the user input, first to input different parameters you need to edit this code, you dont have IO.
	//Second I didnt implement any exception system, this program doesnt have a way to detect semantic erros, specially in the statistics
	//you got to be very careful with your input parameters.

	//IMPORTANT JITTERING SHOULD ONLY BE USED WITH 1 MS BINS AND 5 MS EPOCH. IF YOU WANT TO USE JITTERING WITH OTHER PARAMS YOU NEED TO MODIFY THE SOURCE CODE.

	//Statistician ctor interval is in sec, binsize and epoch is in ms
	bool IsSpontaneous = true; //For Spontanepus activity correlations. if this is false mutlithreading and PREX should be false.
	bool MultiThreading = true; 
	bool PREX = false; //Use first respiration, if its false its gonna use the Stimulus on and off.

	//std::string FileName;

	//std::cout << "File name: ";

	//std::getline(std::cin, FileName);

	//std::string FileName("for_PkC_cs-ss.dat"); //Diff (Cs vs all). 20ms jitter. Block excitation
	//std::string FileName("for_MLI_1.dat"); //Diff (SpMLI vs SpPkC). 3ms jitter. Block excitation
	//std::string FileName("for_MLIsyncOrInh_2or3.dat"); //Same (SpMLI). 3ms jitter. Block inhibition. SYNC
	//std::string FileName("for_MLIsyncOrInh_2or3.dat"); //Same (SpMLI). 3ms jitter. Block excitation. MLI2 inh
	//std::string FileName("for_CS_MLI.dat"); //Diff (Cs vs SpMLI). 10ms jitter. No block. Cs_MLI
	std::string FileName("for_CS_PkC.dat"); //Diff (Cs vs SpMLI). 10ms jitter. No block. Cs_PkC
	//std::string FileName("for_CS_MLI1.dat"); //One2One (Cs vs SpMLI). 10ms jitter. Block excitation. Cs_MLI1
	//std::string FileName("test.dat"); //Name of the File that was created with Matlab code.
	int BinSize = 30; //Samples.1 ms binsize. 1 ms = 30 samples.
	int Epoch = 900; //Samples. Epoch for the analysis. 150 samples = 5ms. 900 samples = 30 ms.
	int JitterResolution = 10;
	uint32_t Interval = 30000; //Seconds. Interval used for statician ctor with PREX enabled.
	uint8_t ResamplingMethod = INTERJITTER; //Select the resampling Method.
	uint8_t StatTest = PERMUTATIONTEST; //Select the Statistics. If PERMUTATIONTEST, ExcZeroLag is ignored.
	int ResampledSets = 1000; //Recommended 100 for shuffle (Burgos-Robles,2017), 1000 for jittering (Fujisawa,2008)
	double ZThreshorPVal = 0.01; //You should calculate this threshold with a two tail Z table. Divide 0.01 / NoBins and then look for the corresponding Z value.
	//p < 0.01 : 3.23 for 8 bins, 3.29 for 10 bins.... p < 0.001 3.84 for 8 bins, 3.89 for 10 bins
	//Put the desire alpha level in this param when using the permutation test.
	bool ExcZeroLag = false; //Important this should only be selected true when the binsize is 1ms if its greater its gonna return garbage.
	
	if (Epoch % BinSize)
	{
		std::cout << " Wrong Epoch or BinSize, not divisible.\n";
		std::cin.get();
		return -1;
	}

	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

	if (PREX)
	{
		std::cout << " Constructing with PREX with " << Interval << " sec interval"<<".\n";
		Statistician SpikeJuggler(FileName, BinSize, Epoch, Interval);
		if (MultiThreading)
		{
			//this is just a warning it doesnt have any effect on the flow of the program.
			std::cout << std::thread::hardware_concurrency() << " concurrent threads are supported.\n"
				<< "You no. of stimuli (Odors) should be less or equal than this number.\n";
			std::cin.get();
			SpikeJuggler.RunThreadPool(ResampledSets,ResamplingMethod, StatTest, ZThreshorPVal, ExcZeroLag);
		}
		else
		{
			std::cout << "Starting in a single thread.\n";
			std::cin.get();
			SpikeJuggler.RunSingleThread(ResampledSets, ResamplingMethod, StatTest, ZThreshorPVal, ExcZeroLag);
		}
	}
	else
	{
		//std::cout << " Constructing with On and Off.\n";

		std::cout << " ResamplingMethod = Interval Jitter\n";
		std::cout << " BinSize = 1ms\n";
		std::cout << " Epoch(Jittering) = �30ms\n";
		std::cout << " Epoch(Permutation test) = �5ms.\n";
		std::cout << " Interval = 3ms.\n\n";




		Statistician SpikeJuggler(FileName, BinSize, Epoch, IsSpontaneous,JitterResolution);
		if (MultiThreading)
		{
			std::cout << std::thread::hardware_concurrency() << " concurrent threads are supported.\n";

			if (std::thread::hardware_concurrency() < 4)
			{
				std::cout << "Multithreading not supported by the system " << "\n";
				std::cin.get();
				return -1;
			}

			std::cout << "Press Enter to start " << "\n";
			std::cin.get();
			SpikeJuggler.RunThreadPool(ResampledSets, ResamplingMethod, StatTest,  ZThreshorPVal, ExcZeroLag);
		}
		else
		{
			std::cout << "Starting in a single thread.\n";
			//std::cin.get();
			SpikeJuggler.RunSingleThread(ResampledSets, ResamplingMethod, StatTest, ZThreshorPVal, ExcZeroLag);
		}
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::chrono::duration<float> duration = end - start;
	

	std::cout << "Finished. Duration in seconds: " << duration.count() << std::endl;
	
	std::cin.get();
}