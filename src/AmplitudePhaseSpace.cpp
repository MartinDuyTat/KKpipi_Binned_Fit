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
  /*phase += KKpipi_Constants::gamma;
  if(phase > TMath::Pi()) {
    phase -= 2*TMath::Pi();
  }
  if(phase < 0.0) {
    return static_cast<int>((phase + TMath::Pi())/(TMath::Pi()/4.0));
  } else if(phase < TMath::Pi()/6) {
    return 4;
  } else if(phase < KKpipi_Constants::gamma) {
    return 5;
  } else if(phase < 2*KKpipi_Constants::gamma - 0.3) {
    return 6;
  } else {
    return 7;
    }*/
  return static_cast<int>((phase + TMath::Pi())/(2*TMath::Pi()/NumberOfBins()));
} 

int AmplitudePhaseSpace::NumberOfBins() const {
  return PhaseSpaceParameterisation::NumberOfBins();
}
