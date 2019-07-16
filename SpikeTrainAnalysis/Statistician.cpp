#include "Statistician.h"

#define REFERENCE 0
#define TARGET 1

Statistician::Statistician(std::string FileName)
	:
	OdorEx(FileName),
	Reference(OdorEx.RDataFile(), OdorEx.GetUnitsRef(), OdorEx.GetRefSizePos(), OdorEx.GetRefTrainPos()),
	Target(OdorEx.RDataFile(), OdorEx.GetUnitsTar(), OdorEx.GetTarSizePos(), OdorEx.GetTarTrainPos()),
	StimLockedSpikesRef(OdorEx.GetStimuli() * OdorEx.GetMagnitudes() * OdorEx.GetTrials() * OdorEx.GetUnitsRef()),
	StimLockedSpikesTar(OdorEx.GetStimuli() * OdorEx.GetMagnitudes() * OdorEx.GetTrials() * OdorEx.GetUnitsTar())
{


}

void Statistician::SetStimLockedSpikes(char RegionType)
{
	if (RegionType == REFERENCE)
	{
		for (auto On = OdorEx.GetStimOn().begin(), Off = OdorEx.GetStimOff().begin(), end = OdorEx.GetStimOn().end();
			On < end;
			++On, ++Off)
		{



			//DataFile->read(reinterpret_cast<char*>(i->data()), i->size() * 8);

		}
	}
	else if (RegionType == TARGET)
	{



	}

}
