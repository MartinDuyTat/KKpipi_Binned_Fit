// Martin Duy Tat 2nd November 2020

#include<algorithm>
#include<vector>
#include<complex>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
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
  int counter = 0;
  std::vector<int> v = binlist.GetEvents(+1);
  while(!std::all_of(v.begin(), v.end(), [&events](int i){return i == events;})) {
    std::vector<TLorentzVector> eventvectors = generator.Generate();
    binlist.AddEvent(Event(eventvectors), +1, events);
    v = binlist.GetEvents(+1);
    counter++;
    if(counter%events == 0) {
      std::cout << "Generated " << counter << " events\n";
    }
  }
  int NumberBins = binlist.NumberBins();
  m_K = std::vector<double>(NumberBins, 0.0);
  m_Kbar = std::vector<double>(NumberBins, 0.0);
  m_c = std::vector<double>(NumberBins, 0.0);
  m_s = std::vector<double>(NumberBins, 0.0);
  Amplitude amplitude("D0toKKpipi.so", "Dbar0toKKpipi.so");
  // Do Monte Carlo integration to find the hadronic parameters
  for(int i = 0; i < NumberBins; i++) {
    std::vector<Event> eventlist_d = binlist.GetBin(i).GetEvents(+1).GetEvents();
    std::vector<Event> eventlist_dbar = binlist.GetBin(i).GetEvents(+1).GetEvents();
    for(int j = 0; j < events; j++) {
      // Calculate amplitude of event
      std::complex<double> amplitude_d = amplitude(eventlist_d[j].GetEvent(), +1);
      std::complex<double> amplitude_dbar = amplitude(eventlist_dbar[j].GetEvent(), -1);
      // Fractional yield in bin i
      m_K[i] += std::norm(amplitude_d);
      m_Kbar[i] += std::norm(amplitude_dbar);
      // Strong Phase difference
      double phase = std::arg(amplitude_d) - std::arg(amplitude_dbar);
      m_c[i] += TMath::Sqrt(std::norm(amplitude_d)*std::norm(amplitude_dbar))*TMath::Cos(phase);
      m_s[i] += TMath::Sqrt(std::norm(amplitude_d)*std::norm(amplitude_dbar))*TMath::Sin(phase);
    }
  }
  double sumK = 0, sumKbar = 0;
  for(int i = 0; i < NumberBins; i++) {
    // Normalise fractional yields so they sum to 1
    sumK += m_K[i];
    sumKbar += m_Kbar[i];
    // Amplitude averaged strong phase variation normalisation
    m_c[i] /= TMath::Sqrt(m_K[i]*m_Kbar[i]);
    m_s[i] /= TMath::Sqrt(m_K[i]*m_Kbar[i]);
  }
  std::transform(m_K.begin(), m_K.end(), m_K.begin(), std::bind(std::divides<double>(), std::placeholders::_1, sumK));
  std::transform(m_Kbar.begin(), m_Kbar.end(), m_Kbar.begin(), std::bind(std::divides<double>(), std::placeholders::_1, sumKbar));
}

DDecayParameters(std::string filename) {
  std::ifstream DDecayFile(filename);
  std::string line;
  std::getline(DDecayFile, line);
  while(std::getline(DDecayFile, line)) {
    std::stringstream ss(line);
    double i, K, Kbar, c, s;
    std::getline(ss, i, ",");
    std::getline(ss, K, ",");
    std::getline(ss, Kbar, ",");
    std::getline(ss, c, ",");
    std::getline(ss, s, "\n");
    m_K.push_back(K);
    m_Kbar.push_back(Kbar);
    m_c.push_back(c);
    m_s.push_back(s);
  }
  DDecayFile.close();
}

void DDecayParameters::saveCSV(filename) const {
  std::ofstream DDecayFile(filename);
  DDecayFile << "i,K_i,Kbar_i,c_i,s_i\n";
  for(int i = 0; i < m_K.size(); i++) {
    DDecayfile << m_K[i] "," << m_Kbar[i] << "," << c_i[i] << "," << s_i[i] << std::endl;
  }
  DDecayFile.close();

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
