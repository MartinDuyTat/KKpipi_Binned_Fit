// Martin Duy Tat 2nd November 2020

#include<vector>
#include<functional>
#include"Bin.h"
#include"BinList.h"
#include"PhaseSpaceParameterisation.h"
#include"Event.h"
#include"TTree.h"
#include"DDecayParameters.h"
#include"CPParameters.h"
#include"TMath.h"

BinList::BinList(PhaseSpaceParameterisation *psp): m_psp(psp), m_bins(std::vector<Bin>(m_psp->NumberOfBins())) {
}

void BinList::AddEvent(Event event, int charge) {
  m_bins[m_psp->WhichBin(event)].AddEvent(event, charge);
}

void BinList::AddEvent(Event event, int charge, int maxevents) {
  int whichbin = m_psp->WhichBin(event);
  if(maxevents > m_bins[whichbin].GetNumberEvents(charge)) {
    m_bins[whichbin].AddEvent(event, charge);
  }
}

void BinList::LoadTTree(TTree *tree, int charge) {
  std::vector<Double_t> four_momentum(16);
  double *p = four_momentum.data();
  tree->SetBranchAddress("_1_K~_Px", p + 0);
  tree->SetBranchAddress("_1_K~_Py", p + 1);
  tree->SetBranchAddress("_1_K~_Pz", p + 2);
  tree->SetBranchAddress("_1_K~_E", p + 3);
  tree->SetBranchAddress("_2_K#_Px", p + 4);
  tree->SetBranchAddress("_2_K#_Py", p + 5);
  tree->SetBranchAddress("_2_K#_Pz", p + 6);
  tree->SetBranchAddress("_2_K#_E", p + 7);
  tree->SetBranchAddress("_3_pi~_Px", p + 8);
  tree->SetBranchAddress("_3_pi~_Py", p + 9);
  tree->SetBranchAddress("_3_pi~_Pz", p + 10);
  tree->SetBranchAddress("_3_pi~_E", p +11);
  tree->SetBranchAddress("_4_pi#_Px", p + 12);
  tree->SetBranchAddress("_4_pi#_Py", p + 13);
  tree->SetBranchAddress("_4_pi#_Pz", p + 14);
  tree->SetBranchAddress("_4_pi#_E", p + 15);
  for(Int_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    this->AddEvent(Event(four_momentum), charge);
  }
}

int BinList::NumberBins() {
  return m_bins.size();
}

std::vector<int> BinList::GetEvents(int charge) const {
  std::vector<int> number_events;
  for(unsigned int i = 0; i < m_bins.size(); i++) {
    number_events.push_back(m_bins[i].GetNumberEvents(charge));
  }
  return number_events;
}

Bin BinList::GetBin(int i) {
  return m_bins[i];
}

void BinList::Predict(const DDecayParameters &ddparameters, const CPParameters &cpparameters, std::vector<double> &BplusEvents, std::vector<double> &BminusEvents, int totalBplus, int totalBminus) {
  double xplus, xminus, yplus, yminus;
  cpparameters.GetCPParameters(xplus, xminus, yplus, yminus);
  std::vector<double> K = ddparameters.GetK();
  std::vector<double> Kbar = ddparameters.GetKbar();
  std::vector<double> c = ddparameters.Getc();
  std::vector<double> s = ddparameters.Gets();
  BplusEvents.resize(m_bins.size());
  BminusEvents.resize(m_bins.size());
  double sumplus = 0, summinus = 0;
  for(unsigned int i = 0; i < m_bins.size(); i++) {
    BplusEvents[i] = Kbar[i] + (xplus*xplus + yplus*yplus)*K[i] + 2*TMath::Sqrt(K[i]*Kbar[i])*(xplus*c[i] - yplus*s[i]);
    BminusEvents[i] = K[i] + (xminus*xminus + yminus*yminus)*Kbar[i] + 2*TMath::Sqrt(K[i]*Kbar[i])*(xminus*c[i] + yminus*s[i]);
    sumplus += BplusEvents[i];
    summinus += BminusEvents[i];
  }
  double normalisationBplus = totalBplus/sumplus;
  double normalisationBminus = totalBminus/summinus;
  std::transform(BplusEvents.begin(), BplusEvents.end(), BplusEvents.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, normalisationBplus));
  std::transform(BminusEvents.begin(), BminusEvents.end(), BminusEvents.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, normalisationBminus));
}
