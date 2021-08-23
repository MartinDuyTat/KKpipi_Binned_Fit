// Martin Duy Tat 12th November 2020

#include"PhaseSpaceParameterisation.h"
#include"Event.h"
#include"TMath.h"
#include"Constants.h"

PhaseSpaceParameterisation::PhaseSpaceParameterisation(const int &bins): m_bins(bins), m_KSVeto(std::pair<double, double>({-1.0, -1.0})) {
}

PhaseSpaceParameterisation::~PhaseSpaceParameterisation() {
}

int PhaseSpaceParameterisation::NumberOfBins() const {
  return m_bins;
}

void PhaseSpaceParameterisation::SetKSVeto(double Lower, double Upper) {
  m_KSVeto = std::pair<double, double>({Lower, Upper});
}

bool PhaseSpaceParameterisation::isKSVeto(const Event &event) const {
  if(m_KSVeto.first <= 0.0 || m_KSVeto.second <= 0.0) {
    return false;
  } else {
    double Mpipi = event.GetInvMass2(2, 3);
    return Mpipi > m_KSVeto.first && Mpipi < m_KSVeto.second;
  }
}
