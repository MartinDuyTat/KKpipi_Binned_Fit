// Martin Duy Tat 9th December 2020

#include<vector>
#include<complex>
#include"PhaseSpaceParameterisation.h"
#include"AmplitudePhaseSpace.h"
#include"Amplitude.h"
#include"Event.h"

AmplitudePhaseSpace::AmplitudePhaseSpace(int bins): PhaseSpaceParameterisation(bins) {
}

int AmplitudePhaseSpace::WhichBin(const Event &event) const {
  std::vector<double> EventVector = event.GetEventVector();
  double phase = std::arg(m_amplitude(EventVector, +1)*std::conj(m_amplitude(EventVector, -1)));
  return static_cast<int>((phase + TMath::Pi())/(2*TMath::Pi()/NumberOfBins()));
} 

int AmplitudePhaseSpace::NumberOfBins() const {
  return PhaseSpaceParameterisation::NumberOfBins();
}
