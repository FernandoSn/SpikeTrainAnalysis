#include <iostream>
#include "Statistician.h"
#include <algorithm>
#include <chrono>
#include <random>


int main()
{
	//Statistician ctor interval is in sec, binsize and epoch is in ms


	std::vector<int> v{ 3, 1, 4, 1, 5, 9, 2, 6 };
	auto bounds = std::minmax_element(v.begin(), v.end());

	std::random_device r;
	std::default_random_engine generator(r());
	std::uniform_int_distribution<int> distribution(0, 0);
	int asd = distribution(generator);
	double rrrr = 99999999.2;
	auto ty = std::for_each(v.begin(), v.begin(), [](auto& n) { n *= 2; });
	//auto ff = *ty;
	sizeof(distribution);


	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
	Statistician SpikeJuggler("Prueba.dat", 1.0,25,100);
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

	std::cout << duration.count() << std::endl;
	
	std::cin.get();


}