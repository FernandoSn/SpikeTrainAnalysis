#include <iostream>
#include "BrainRegion.h"

#define REFERENCE 6
#define TARGET 8


int main()
{


	std::cout << "hola" << std::endl;

	BrainRegion AON("Prueba.dat",REFERENCE);
	BrainRegion PfCx("Prueba.dat", TARGET);



	std::cin.get();

}