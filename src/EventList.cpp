// Martin Duy Tat 2nd November 2020

#include<vector>
#include<string>
#include"TFile.h"
#include"TTree.h"
#include"EventList.h"
#include"Event.h"
#include"KKpipiMath.h"

EventList::EventList() {
}

void EventList::AddEvent(const Event &event) {
  m_eventlist.push_back(event);
}

int EventList::NumberEvents() const {
  return m_eventlist.size();
}

const std::vector<Event>& EventList::GetEvents() const {
  return m_eventlist;
}

const Event& EventList::GetEvent(const int &i) const {
  return m_eventlist[i];
}

void EventList::LoadTree(const std::string &filename, const std::string &filename_weight) {
  TFile f(filename.c_str(), "READ");
  TTree *tree = nullptr;
  f.GetObject("FlatEvents", tree);
  std::vector<double> X(5);
  std::complex<double> *D_amplitude = nullptr, *Dbar_amplitude = nullptr;
  double weight = 1.0;
  tree->SetBranchAddress("x1", X.data() + 0);
  tree->SetBranchAddress("x2", X.data() + 1);
  tree->SetBranchAddress("x3", X.data() + 2);
  tree->SetBranchAddress("x4", X.data() + 3);
  tree->SetBranchAddress("x5", X.data() + 4);
  tree->SetBranchAddress("D_amplitude", &D_amplitude);
  tree->SetBranchAddress("Dbar_amplitude", &Dbar_amplitude);
  if(filename_weight != "") {
    tree->AddFriend("FlatEventsWeights", filename_weight.c_str());
    tree->SetBranchAddress("weight", &weight);
  }
  m_eventlist.reserve(tree->GetEntries());
  const auto N = tree->GetEntries();
  for(int i = 0; i < N; i++) {
    tree->GetEntry(i);
    const std::vector<double> momenta = KKpipiMath::ConvertXToMomenta(X);
    this->AddEvent(Event(momenta, *D_amplitude, *Dbar_amplitude, weight));
  }
}

void EventList::Clear() {
  std::vector<Event>().swap(m_eventlist);
}
