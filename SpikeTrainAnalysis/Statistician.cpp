#include "Statistician.h"

Statistician::Statistician(std::string FileName)
	:
	OdorEx(FileName),
	Reference(OdorEx.RDataFile(), OdorEx.GetUnitsRef(), OdorEx.GetRefSizePos(), OdorEx.GetRefTrainPos()),
	Target(OdorEx.RDataFile(), OdorEx.GetUnitsTar(), OdorEx.GetTarSizePos(), OdorEx.GetTarTrainPos())
{


}
