#include <iostream>
#include "Statistician.h"
#include <algorithm>
#include <chrono>


int main()
{
	int a = 25;
	double c = 1000.0;
	double b = 25 / 1000;


	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
	Statistician SpikeJuggler("Prueba.dat", 1.0);
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