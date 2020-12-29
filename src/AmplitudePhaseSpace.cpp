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
  std::complex<double> D_amplitude, Dbar_amplitude;
  double phase, rD;
  if(m_ReadAmplitudeFromEvent) {    
    event.GetAmplitudes(D_amplitude, Dbar_amplitude);
  } else {
    D_amplitude = m_amplitude(EventVector, +1);
    Dbar_amplitude = m_amplitude(EventVector, -1);
  }
  phase = std::arg(D_amplitude*std::conj(Dbar_amplitude));
  rD = std::norm(D_amplitude/Dbar_amplitude);
  int BinNumber = static_cast<int>((phase + TMath::Pi())/(2*TMath::Pi()/NumberOfBins())) + 1;
  if(rD < 1) {
    return BinNumber;
  } else {
    return -(NumberOfBins() + 1 - BinNumber);
  }
} 

int AmplitudePhaseSpace::NumberOfBins() const {
  return PhaseSpaceParameterisation::NumberOfBins();
}
