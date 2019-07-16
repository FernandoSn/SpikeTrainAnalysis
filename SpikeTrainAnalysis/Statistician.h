#pragma once

#include "Experiment.h"
#include "BrainRegion.h"

class Statistician
{
public:

	Statistician(std::string FileName);

private:

	Experiment OdorEx;
	BrainRegion Reference;
	BrainRegion Target;

	

};