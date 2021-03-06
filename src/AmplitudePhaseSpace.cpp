// Martin Duy Tat 9th December 2020

#include<vector>
#include<complex>
#include"TMath.h"
#include"PhaseSpaceParameterisation.h"
#include"AmplitudePhaseSpace.h"
#include"Amplitude.h"
#include"Event.h"
#include"Constants.h"

AmplitudePhaseSpace::AmplitudePhaseSpace(const int &bins, bool rDBinning): PhaseSpaceParameterisation(bins), m_rDBinning(rDBinning) {
}

AmplitudePhaseSpace::~AmplitudePhaseSpace() {
}

void AmplitudePhaseSpace::ReadAmplitudeFromEvent(bool TrueIfReadFromEvent) {
  m_ReadAmplitudeFromEvent = TrueIfReadFromEvent;
}

void AmplitudePhaseSpace::UseVariableBinWidths(bool VariableBinWidths) {
  if(VariableBinWidths && static_cast<int>(m_BinMap.size()) == NumberOfBins()) {
    m_UseVariableBinWidths = true;
  } else {
    m_UseVariableBinWidths = false;
  }
}

void AmplitudePhaseSpace::SetBinEdges(const std::vector<double> &BinEdges) {
  m_BinMap.clear();
  int BinNumberCounter = 1;
  for(auto it = BinEdges.rbegin(); it != BinEdges.rend(); it++) {
    m_BinMap.insert({-*it, BinNumberCounter});
    ++BinNumberCounter;
  }
  m_BinMap.insert({0.0, BinNumberCounter});
  ++BinNumberCounter;
  for(const auto &edge : BinEdges) {
    m_BinMap.insert({edge, BinNumberCounter});
    ++BinNumberCounter;
  }
  m_BinMap.insert({TMath::Pi(), BinNumberCounter});
}

int AmplitudePhaseSpace::WhichBin(const Event &event) const {
  if(isKSVeto(event)) {
    return 0;
  }
  std::vector<double> EventVector = event.GetEventVector();
  std::complex<double> D_amplitude, Dbar_amplitude;
  double phase, rD;
  if(m_ReadAmplitudeFromEvent) {    
    event.GetAmplitudes(D_amplitude, Dbar_amplitude);
  } else {
    D_amplitude = m_amplitude(EventVector, +1);
    Dbar_amplitude = m_amplitude(EventVector, -1);
  }
  if(TMath::IsNaN(D_amplitude.real()) || TMath::IsNaN(D_amplitude.imag()) || TMath::IsNaN(Dbar_amplitude.real()) || TMath::IsNaN(Dbar_amplitude.imag())) {
    return 0;
  }
  phase = std::arg(D_amplitude*std::conj(Dbar_amplitude));
  rD = std::norm(D_amplitude/Dbar_amplitude);
  int BinNumber;
  if(m_UseVariableBinWidths) {
    BinNumber = m_BinMap.lower_bound(phase)->second;
  } else {
    BinNumber = static_cast<int>((phase + TMath::Pi())/(2*TMath::Pi()/NumberOfBins())) + 1;
  }
  if(!m_rDBinning) {
    if(rD < 1) {
      return BinNumber;
    } else {
      return -(NumberOfBins() + 1 - BinNumber);
    }
  } else {
    double ln_rD = 0.5*TMath::Log(rD);
    if(ln_rD < 0.0 && ln_rD > -0.5) {
      return BinNumber;
    } else if(ln_rD < -0.5) {
      return BinNumber + NumberOfBins();
    } else if(ln_rD > 0.0 && ln_rD < 0.5) {
      return -(NumberOfBins() + 1 - BinNumber);
    } else {
      return -(2*NumberOfBins() + 1 - BinNumber);
    }
  }
} 

int AmplitudePhaseSpace::NumberOfBins() const {
  if(m_rDBinning) {
    return PhaseSpaceParameterisation::NumberOfBins()/2;
  } else {
    return PhaseSpaceParameterisation::NumberOfBins();
  }
}
