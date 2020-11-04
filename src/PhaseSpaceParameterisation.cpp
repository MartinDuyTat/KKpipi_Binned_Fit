// Martin Duy Tat 2nd November 2020

#include<vector>
#include"PhaseSpaceParameterisation.h"
#include"Event.h"

PhaseSpaceParameterisation::PhaseSpaceParameterisation(): m_bins(4) {
}

int PhaseSpaceParameterisation::WhichBin(const Event &event) {
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

int PhaseSpaceParameterisation::NumberOfBins() {
  return m_bins;
}
