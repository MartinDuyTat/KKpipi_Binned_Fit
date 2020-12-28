// Martin Duy Tat 9th December 2020

#include<vector>
#include<complex>
#include"TMath.h"
#include"PhaseSpaceParameterisation.h"
#include"AmplitudePhaseSpace.h"
#include"Amplitude.h"
#include"Event.h"
#include"Constants.h"

AmplitudePhaseSpace::AmplitudePhaseSpace(const int &bins): PhaseSpaceParameterisation(bins) {
}

AmplitudePhaseSpace::~AmplitudePhaseSpace() {
}

void AmplitudePhaseSpace::ReadAmplitudeFromEvent(bool TrueIfReadFromEvent) {
  m_ReadAmplitudeFromEvent = TrueIfReadFromEvent;
}

int AmplitudePhaseSpace::WhichBin(const Event &event) const {
  std::vector<double> EventVector = event.GetEventVector();
  double phase;
  if(m_ReadAmplitudeFromEvent) {
    std::complex<double> D_amplitude, Dbar_amplitude;
    event.GetAmplitudes(D_amplitude, Dbar_amplitude);
    phase = std::arg(D_amplitude*std::conj(Dbar_amplitude));
  } else {
    phase = std::arg(m_amplitude(EventVector, +1)*std::conj(m_amplitude(EventVector, -1)));
  }
  if(phase > 0) {
    return static_cast<int>(phase/(TMath::Pi()/NumberOfBins())) + 1;
  } else {
    return -(static_cast<int>(-phase/(TMath::Pi()/NumberOfBins())) + 1);
  }
} 

int AmplitudePhaseSpace::NumberOfBins() const {
  return PhaseSpaceParameterisation::NumberOfBins();
}
