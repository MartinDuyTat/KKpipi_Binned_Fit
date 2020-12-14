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

void BinList::AddEvent(const Event &event, const int &charge) {
  m_bins[m_psp->WhichBin(event)].AddEvent(event, charge);
}

void BinList::AddEvent(const Event &event, const int &charge, const int &maxevents) {
  int whichbin = m_psp->WhichBin(event);
  if(maxevents > m_bins[whichbin].GetNumberEvents(charge)) {
    m_bins[whichbin].AddEvent(event, charge);
  }
}

void BinList::LoadTTree(TTree *tree, const int &charge, const int &StartEvent, const int &TotalEvents) {
  std::vector<Double_t> four_momentum(16);
  double *p = four_momentum.data();
  for(int i = 0; i < 16; i++) {
    std::string address = tree->GetListOfBranches()[0][i]->GetName();
    if(i%4 == 0) {
      tree->SetBranchAddress(address.c_str(), p + i + 3);
    } else {
      tree->SetBranchAddress(address.c_str(), p + i - 1);
    }
  }
  for(Int_t i = StartEvent; i < ((TotalEvents  == -1) ? tree->GetEntries() : (StartEvent + TotalEvents)); i++) {
    tree->GetEntry(i);
    this->AddEvent(Event(four_momentum), charge);
  }
}

int BinList::NumberBins() {
  return m_bins.size();
}

std::vector<int> BinList::GetEvents(const int &charge) const {
  std::vector<int> number_events;
  for(unsigned int i = 0; i < m_bins.size(); i++) {
    number_events.push_back(m_bins[i].GetNumberEvents(charge));
  }
  return number_events;
}

Bin BinList::GetBin(const int &i) {
  return m_bins[i];
}

void BinList::Predict(const DDecayParameters &ddparameters, const CPParameters &cpparameters, std::vector<double> &BplusEvents, std::vector<double> &BminusEvents, const int &totalBplus, const int &totalBminus) {
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
