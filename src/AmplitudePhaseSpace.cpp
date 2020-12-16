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

int AmplitudePhaseSpace::WhichBin(const Event &event) const {
  std::vector<double> EventVector = event.GetEventVector();
  double phase = std::arg(m_amplitude(EventVector, +1)*std::conj(m_amplitude(EventVector, -1)));
  phase -= KKpipi_Constants::gamma - TMath::Pi()/6;
  if(phase < -TMath::Pi()) {
    phase += 2*TMath::Pi();
  }
  if(phase < -7*TMath::Pi()/9) {
    return 0;
  } else if(phase >= -7*TMath::Pi()/9 && phase < -5*TMath::Pi()/9) {
    return 1;
  } else if(phase >= -5*TMath::Pi()/9 && phase < -TMath::Pi()/3) {
    return 2;
  } else if(phase >= -TMath::Pi()/3 && phase < -TMath::Pi()/6) {
    return 3;
  } else if(phase >= -TMath::Pi()/6 && phase < 0.0) {
    return 4;
  } else if(phase >= 0.0 && phase < TMath::Pi()/3) {
    return 5;
  } else if(phase >= TMath::Pi()/3 && phase < 2*TMath::Pi()/3) {
    return 6;
  } else {
    return 7;
  }
  //return static_cast<int>((phase + TMath::Pi())/(2*TMath::Pi()/NumberOfBins()));
} 

int AmplitudePhaseSpace::NumberOfBins() const {
  return PhaseSpaceParameterisation::NumberOfBins();
}
