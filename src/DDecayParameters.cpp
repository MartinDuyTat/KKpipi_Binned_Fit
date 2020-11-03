// Martin Duy Tat 2nd November 2020

#include<algorithm>
#include<vector>
#include<complex>
#include"DDecayParameters.h"
#include"BinList.h"
#include"Generator.h"
#include"TLorentzVector.h"
#include"Event.h"
#include"EventList.h"
#include"Bin.h"
#include"TMath.h"
#include"Amplitude.h"

DDecayParameters::DDecayParameters(const PhaseSpaceParameterisation &psp, const double &mass_parent, const Double_t *mass_decay, int events) {
  Generator generator(mass_parent, mass_decay, 4);
  BinList binlist(psp);
  // Generate events until all bins have "events" number of events
  // Use different events for K/Kbar and c/s
  while(std::all_of(binlist.GetEvents(+1).begin(), binlist.GetEvents(-1).end(), [&events](int i){return i == events;})) {
    std::vector<TLorentzVector> eventvectors = generator.Generate();
    binlist.AddEvent(Event(eventvectors), +1, events);
  }		    
  while(std::all_of(binlist.GetEvents(-1).begin(), binlist.GetEvents(-1).end(), [&events](int i){return i == events;})) {
    std::vector<TLorentzVector> eventvectors = generator.Generate();
    binlist.AddEvent(Event(eventvectors), -1, events);
  }
  int NumberBins = binlist.NumberBins();
  m_K = std::vector<double>(NumberBins, 0.0);
  m_Kbar = std::vector<double>(NumberBins, 0.0);
  m_c = std::vector<double>(NumberBins, 0.0);
  m_s = std::vector<double>(NumberBins, 0.0);
  Amplitude amplitude("D0toKKpipi.so", "Dbar0toKKpipi.so");
  for(int i = 0; i < NumberBins; i++) {
    std::vector<Event> eventlist_d = binlist.GetBin(i).GetEvents(-1).GetEvents();
    std::vector<Event> eventlist_dbar = binlist.GetBin(i).GetEvents(+1).GetEvents();
    for(int j = 0; j < events; j++) {
      std::complex<double> amplitude_d = amplitude(eventlist_d[j].GetEvent(), +1);
      std::complex<double> amplitude_dbar = amplitude(eventlist_dbar[j].GetEvent(), +1);
      m_K[i] += std::norm(amplitude_d)/events;
      m_Kbar[i] += std::norm(amplitude_dbar)/events;
      std::complex<double>amplitude_ddbar = amplitude_d/amplitude_dbar;
      m_c[i] += TMath::Sqrt(std::norm(amplitude_d))*TMath::Sqrt(std::norm(amplitude_dbar))*amplitude_ddbar.real()/events;
      m_s[i] += TMath::Sqrt(std::norm(amplitude_d))*TMath::Sqrt(std::norm(amplitude_dbar))*amplitude_ddbar.imag()/events;
    }
  }
  double sumK = 0, sumKbar = 0;
  for(int i = 0; i < NumberBins; i++) {
    sumK += m_K[i];
    sumKbar += m_Kbar[i];
    m_c[i] /= TMath::Sqrt(m_K[i]*m_Kbar[i]);
    m_s[i] /= TMath::Sqrt(m_K[i]*m_Kbar[i]);
  }
  std::transform(m_K.begin(), m_K.end(), m_K.begin(), std::bind(std::divides<double>(), std::placeholders::_1, sumK));
  std::transform(m_Kbar.begin(), m_Kbar.end(), m_Kbar.begin(), std::bind(std::divides<double>(), std::placeholders::_1, sumKbar));
}

std::vector<double> DDecayParameters::GetK() const {
  return m_K;
}
std::vector<double> DDecayParameters::GetKbar() const {
  return m_Kbar;
}

std::vector<double> DDecayParameters::Getc() const {
  return m_c;
}

std::vector<double> DDecayParameters::Gets() const {
  return m_s;
}
