// Martin Duy Tat 2nd November 2020

#include<vector>
#include"Bin.h"
#include"BinList.h"
#include"PhaseSpaceParameterisation.h"
#include"Event.h"
#include"TTree.h"
#include"TMath.h"

BinList::BinList(PhaseSpaceParameterisation *psp, bool SaveEvents): m_psp(psp), m_bins(std::vector<Bin>(m_psp->NumberOfBins())), m_CPbins(std::vector<Bin>(m_psp->NumberOfBins())), m_SaveEvents(SaveEvents), m_BplusEvents(std::vector<int>(m_psp->NumberOfBins())), m_BminusEvents(std::vector<int>(m_psp->NumberOfBins())), m_BplusCPEvents(std::vector<int>(m_psp->NumberOfBins())), m_BminusCPEvents(std::vector<int>(m_psp->NumberOfBins())) {
}

void BinList::AddEvent(const Event &event, const int &charge) {
  int BinNumber = m_psp->WhichBin(event);
  if(BinNumber == 0) {
    return;
  }
  if(charge == +1) {
    if(BinNumber > 0) {
      m_BplusEvents[BinNumber - 1] += 1;
    } else {
      m_BplusCPEvents[-BinNumber - 1] += 1;
    }
  } else if (charge == -1) {
    if(BinNumber > 0) {
      m_BminusEvents[BinNumber - 1] += 1;
    } else {
      m_BminusCPEvents[-BinNumber - 1] += 1;
    }
  }
  if(m_SaveEvents) {
    if(BinNumber > 0) {
      m_bins[BinNumber - 1].AddEvent(event, charge);
    } else {
      m_CPbins[-BinNumber - 1].AddEvent(event, charge);
    }
  }
}

void BinList::AddEvent(const Event &event, const int &charge, const int &maxevents) {
  int BinNumber = m_psp->WhichBin(event);
  if(BinNumber == 0) {
    return;
  }
  if(charge == +1) {
    if(BinNumber > 0) {
      m_BplusEvents[BinNumber - 1] += 1;
    } else {
      m_BplusCPEvents[-BinNumber - 1] += 1;
    }
  } else if (charge == -1) {
    if(BinNumber > 0) {
      m_BminusEvents[BinNumber - 1] += 1;
    } else {
      m_BminusCPEvents[-BinNumber - 1] += 1;
    }
  }
  if(m_SaveEvents && maxevents > m_bins[BinNumber].GetNumberEvents(charge)) {
    if(BinNumber > 0) {
      m_bins[BinNumber - 1].AddEvent(event, charge);
    } else {
      m_CPbins[-BinNumber - 1].AddEvent(event, charge);
    }
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

std::vector<int> BinList::GetEvents(const int &charge, const int &cp) const {
  if(charge == +1) {
    if(cp == +1) {
      return m_BplusEvents;
    } else if (cp == -1) {
      return m_BplusCPEvents;
    }
  } else if(charge == -1) {
    if(cp == +1) {
      return m_BminusEvents;
    } else if(cp == -1) {
      return m_BminusCPEvents;
    }
  }
  return {};
}

Bin BinList::GetBin(const int &i, const int &cp) {
  if(cp == +1) {
    return m_bins[i];
  } else if(cp == -1) {
    return m_CPbins[i];
  } else {
    return Bin();
  }
}
