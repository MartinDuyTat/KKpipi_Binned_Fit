// Martin Duy Tat 12th November 2020

#include"RectangularPhaseSpace.h"
#include"PhaseSpaceParameterisation.h"

RectangularPhaseSpace::RectangularPhaseSpace(): PhaseSpaceParameterisation(1) {
}

int PhaseSpaceParameterisation::WhichBin(const Event &event) {
  if(event.GetInvMass2(0, 1) > 0) {
    return 0;
  } else {
    return 0;
  }
}

int RectangularPhaseSpace::NumberOfBins() const {
  return PhaseSpaceParameterisation::NumberOfBins();
}
