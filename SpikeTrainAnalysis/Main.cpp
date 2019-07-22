#include <iostream>
#include <chrono>
#include <thread>
#include <string>
#include <fstream>

#include "Statistician.h"

int main()
{
	//This is not professional coding. to my taste class are well designed maybe the thread pool or MasterCorr on Statistician
	//needs more work because a have a ton of code in one method but its acceptable.
	//The part that is really bad is the user input, first to input different parameters you need to edit this code, you dont have IO.
	//Second I didnt implement any exception system, this program doesnt have a way to detect semantic erros, specially in the statistics
	//you got to be very careful with your input parameters.

	//Statistician ctor interval is in sec, binsize and epoch is in ms

	bool MultiThreading = false;
	std::string FileName("Prueba.dat"); //Name of the File that was created with Matlab code.
	int BinSize = 25; //Miliseconds.
	int Epoch = 100; //Miliseconds. Epoch for the analysis.
	bool PREX = true; //Use first respiration, if its false its gonna use the Stimulus on and off.
	double Interval = 1.0; //Seconds. Interval used for statician ctor with PREX enabled.
	unsigned char ResamplingMethod = SHUFFLING; //Select the resampling Method.
	int ResampledSets = 100; //Recommended 100 for shuffle (Burgos-Robles,2017), 1000 for jittering (Fujisawa,2008)
	double ZThresh = 3.23; //You should calculate this threshold with a two tail Z table. Divide 0.01 / NoBins and then look for the corresponding Z value.
	//3.23 for 8 bins, 3.29 for 10 bins
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
		Statistician SpikeJuggler(FileName, Interval, BinSize, Epoch);
		if (MultiThreading)
		{
			//this is just a warning it doesnt have any effect on the flow of the program.
			std::cout << std::thread::hardware_concurrency() << " concurrent threads are supported.\n"
				<< "You no. of stimuli (Odors) should be less or equal than this number.\n";
			std::cin.get();
			SpikeJuggler.RunThreadPool(ResampledSets,ResamplingMethod,ZThresh,ExcZeroLag);
		}
		else
		{
			std::cout << "Starting in a single thread.\n";
			std::cin.get();
			SpikeJuggler.MasterSpikeCrossCorr(ResampledSets, ResamplingMethod, ZThresh, ExcZeroLag);
		}
	}
	else
	{
		std::cout << " Constructing with On and Off.\n";
		Statistician SpikeJuggler(FileName, BinSize, Epoch);
		if (MultiThreading)
		{
			std::cout << std::thread::hardware_concurrency() << " concurrent threads are supported.\n"
				<< "You no. of stimuli (Odors) should be less or equal than this number\n";
			std::cin.get();
			SpikeJuggler.RunThreadPool(ResampledSets, ResamplingMethod, ZThresh, ExcZeroLag);
		}
		else
		{
			std::cout << "Starting in a single thread.\n";
			std::cin.get();
			SpikeJuggler.MasterSpikeCrossCorr(ResampledSets, ResamplingMethod, ZThresh, ExcZeroLag);
		}
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::chrono::duration<float> duration = end - start;
	

	std::cout << "Finished. Duration in seconds: " << duration.count() << std::endl;
	
	std::cin.get();
}