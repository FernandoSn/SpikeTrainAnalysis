#include <iostream>
#include "Statistician.h"
#include <algorithm>
#include <chrono>
#include <random>
#include <thread>
#include <string>
#include <fstream>


int main()
{
	//This is not professional coding. to my taste class are well designed maybe the thread pool or MasterCorr on Statistician
	//needs more work because a have a ton of code in one method but its acceptable.
	//The part that is really bad is the user input, first to input different parameters you need to edit this code, you dont have IO.
	//Second I didnt implement any exception system, this program doesnt have a way to detect semantic erros, specially in the statistics
	//you got to be very careful with your input parameters.


	//std::vector<unsigned int> asd(8,12);
	//Statistician ctor interval is in sec, binsize and epoch is in ms

	unsigned int n = std::thread::hardware_concurrency();
	std::cout << n << " concurrent threads are supported.\n";

	std::vector<int> v2{ 3, 1, 4, 1, 5, 9, 2, 6 };
	std::vector<int> v{1, 1, 1, 1, 1, 1, 1, 1};

	v.resize(v2.size()+1);


	//auto ty = std::for_each(v2.begin(), v2.end(), [](auto& n) { n *= 2; });
	std::transform(v2.begin(), v2.end(), v.begin(), v2.begin(),[](auto& n, auto& m) { return n + m; });

	std::transform(v2.begin(), v2.end(), v.begin(), v2.begin(), [](auto& n, auto& m) { return 1; });

	v.erase(std::copy_if(v2.begin(),
		v2.end(),
		v.begin(),
		[](double SpikeTime) {return SpikeTime == 9; }) //Lambda as predicate for the algorithm.
		, v.end());

	v.shrink_to_fit();


	auto bounds = std::minmax_element(v.begin(), v.end());

	std::random_device r;
	std::default_random_engine generator(r());
	std::uniform_int_distribution<int> distribution(0, 0);
	int asd = distribution(generator);
	double rrrr = 99999999.2;
	//auto ff = *ty;
	sizeof(distribution);


	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
	//Statistician SpikeJuggler("Prueba.dat", 1.0,25,100);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::chrono::duration<float> duration = end - start;


	/*std::vector<int> Unit{1,1,1,1,6,57,68,6,1,1,1,56,7,1,77};
	std::vector<int> dest;
	dest.resize(10);
	int thresh = 2;
	int thresh2 = 69;


	dest.erase(std::copy_if(Unit.begin(),
		Unit.end(),
		dest.begin(),
		[thresh,thresh2](int u) {return u > thresh && u < thresh2; }),dest.end());

	dest.shrink_to_fit();*/

	int figl = 10;
	int figl2 = -figl;

	bool rtt = 2 < figl;

	std::cout << duration.count() << std::endl;
	
	std::cin.get();


}