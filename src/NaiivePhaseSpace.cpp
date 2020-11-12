// Martin Duy Tat 2nd November 2020

#include<vector>
#include"NaiivePhaseSpace.h"
#include"Event.h"
#include"PhaseSpaceParameterisation.h"

NaiivePhaseSpace::NaiivePhaseSpace(): PhaseSpaceParameterisation(4) {
}

int NaiivePhaseSpace::WhichBin(const Event &event) {
  std::vector<double> momenta = event.GetEvent();
  if(momenta[3] > momenta[7] && momenta[11] > momenta[15]) {
    return 0;
  } else if(momenta[3] > momenta[7] && momenta[11] < momenta[15]) {
    return 1;
  } else if(momenta[3] < momenta[7] && momenta[11] > momenta[15]) {
    return 2;
  } else {
    return 3;
  }
}

int NaiivePhaseSpace::NumberOfBins() const {
  return PhaseSpaceParameterisation::NumberOfBins();
}
