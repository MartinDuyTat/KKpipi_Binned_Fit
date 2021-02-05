// Martin Duy Tat 12th November 2020

#include"PhaseSpaceParameterisation.h"
#include"Event.h"
#include"TMath.h"
#include"Constants.h"

PhaseSpaceParameterisation::PhaseSpaceParameterisation(const int &bins): m_bins(bins) {
}

PhaseSpaceParameterisation::~PhaseSpaceParameterisation() {
}

int PhaseSpaceParameterisation::NumberOfBins() const {
  return m_bins;
}

void PhaseSpaceParameterisation::SetKSVeto(double veto) {
  m_KSVeto = veto;
}

bool PhaseSpaceParameterisation::isKSVeto(const Event &event) const {
  if(m_KSVeto <= 0.0) {
    return false;
  } else {
    return TMath::Abs(event.GetInvMass2(2, 3) - KKpipi_Constants::MASS_KS) < m_KSVeto;
  }
}
