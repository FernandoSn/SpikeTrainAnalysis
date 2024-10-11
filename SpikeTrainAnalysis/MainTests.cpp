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
	//std::vector<int> from_vector(10);
	//std::iota(from_vector.begin(), from_vector.end(), 0);

	//for (auto it = from_vector.begin(); it < from_vector.end(); ++it)
	//{
	//	std::cout << *it << "\n";

	//}
	//std::cout << (from_vector.end()+1<from_vector.end()) << "\n";

	std::vector<int> from_vector = {12 ,26 ,123, 56 ,87 ,2354, 877};

	auto itb = from_vector.begin();
	auto ite = from_vector.end();

	int LeadCount = std::accumulate(itb, ite, 0);
	auto asd = std::find(itb, ite, 1);

	std::cout << LeadCount << "\n";

	std::cout << *itb << "\n";
	std::cout << *(ite-1) << "\n";



	std::cin.get();
}