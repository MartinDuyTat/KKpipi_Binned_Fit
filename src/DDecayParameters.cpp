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
#include"TCanvas.h"
#include"TGraph.h"
#include"TAxis.h"

DDecayParameters::DDecayParameters(PhaseSpaceParameterisation *psp, const double &mass_parent, const Double_t *mass_decay, int events) {
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

DDecayParameters::DDecayParameters(std::string filename) {
  std::ifstream DDecayFile(filename);
  std::string line;
  std::getline(DDecayFile, line);
  while(std::getline(DDecayFile, line)) {
    std::stringstream ss(line);
    int i;
    double K, Kbar, c, s;
    ss >> i;
    ss.ignore();
    ss >> K;
    ss.ignore();
    ss >> Kbar;
    ss.ignore();
    ss >> c;
    ss.ignore();
    ss >> s;
    m_K.push_back(K);
    m_Kbar.push_back(Kbar);
    m_c.push_back(c);
    m_s.push_back(s);
  }
  DDecayFile.close();
}

void DDecayParameters::SaveCSV(std::string filename) const {
  std::ofstream DDecayFile(filename);
  DDecayFile << "i,K_i,Kbar_i,c_i,s_i\n";
  for(unsigned int i = 0; i < m_K.size(); i++) {
    DDecayFile << i << "," <<  m_K[i] << "," << m_Kbar[i] << "," << m_c[i] << "," << m_s[i] << std::endl;
  }
  DDecayFile.close();
}

void DDecayParameters::PlotParameters(std::string filename_cs, std::string filename_K) {
  std::vector<Double_t> circle_x(101), circle_y(101);
  for(int i = 0; i < 101; i++) {
    circle_x[i] = TMath::Cos(2*TMath::Pi()*i/100);
    circle_y[i] = TMath::Sin(2*TMath::Pi()*i/100);
  }
  TGraph *gr1 = new TGraph(this->Getc().size(), this->Gets().data(), this->Getc().data());
  TGraph *circle = new TGraph(101, circle_x.data(), circle_y.data());
  TCanvas *c1 = new TCanvas("s_vs_c", "s_i vs c_i", 700, 700);
  circle->Draw("AL");
  gr1->Draw("P");
  circle->GetXaxis()->SetTitle("s_{i}");
  circle->GetYaxis()->SetTitle("c_{i}");
  circle->GetXaxis()->SetRangeUser(-1.0, 1.0);
  circle->GetYaxis()->SetRangeUser(-1.0, 1.0);
  circle->SetTitle("Plot of c_i vs s_i");
  gr1->SetMarkerStyle(kFullDotLarge);
  c1->SetLeftMargin(0.11);
  circle->Draw("AL");
  gr1->Draw("P");
  c1->Update();
  c1->SaveAs(filename_cs.c_str());
  std::vector<double> binning(this->GetK().size());
  for(unsigned int i = 0; i < binning.size(); i++) {
    binning[i] = i;
  }
  TGraph *gr2 = new TGraph(binning.size(), binning.data(), this->GetK().data());
  TCanvas *c2 = new TCanvas("K", "K_i and Kbar_i");
  gr2->Draw("*A");
  gr2->GetXaxis()->SetTitle("Bin number");
  gr2->GetYaxis()->SetTitle("Fractional yield");
  gr2->GetXaxis()->SetLimits(-0.5, 3.5);
  gr2->GetYaxis()->SetRangeUser(0.0, 0.4);
  gr2->SetTitle("Fractional yields");
  gr2->SetMarkerStyle(kFullDotLarge);
  gr2->Draw("AP");
  c2->Update();
  c2->SaveAs(filename_K.c_str());
  delete gr1;
  delete gr2;
  delete circle;
  delete c1;
  delete c2;
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
