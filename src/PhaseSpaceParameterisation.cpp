// Martin Duy Tat 12th November 2020

#include"PhaseSpaceParameterisation.h"

PhaseSpaceParameterisation::PhaseSpaceParameterisation(int bins): m_bins(bins) {
}

PhaseSpaceParameterisation::~PhaseSpaceParameterisation() {
}

int PhaseSpaceParameterisation::NumberOfBins() const {
  return m_bins;
}
