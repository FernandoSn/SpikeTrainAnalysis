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
	std::vector<int> from_vector(10);
	//std::iota(from_vector.begin(), from_vector.end(), 0);

	for (auto it = from_vector.begin(); it < from_vector.end(); ++it)
	{
		std::cout << *it << "\n";

	}
	std::cin.get();
}