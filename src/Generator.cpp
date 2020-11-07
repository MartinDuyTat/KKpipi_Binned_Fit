// Martin Duy Tat 2nd November

#include<vector>
#include"Generator.h"
#include"TGenPhaseSpace.h"
#include"TLorentzVector.h"
#include"TRandom3.h"

Generator::Generator(const Double_t &mass_parent, const Double_t *mass_decay, Int_t particles) {
  m_phasespace = TGenPhaseSpace();
  TLorentzVector P(0.0, 0.0, 0.0, mass_parent);
  m_phasespace.SetDecay(P, particles, mass_decay);
}

std::vector<TLorentzVector> Generator::Generate() {
  Double_t weight = 1.0;
  TRandom3 random_generator(0);
  do {
    weight = m_phasespace.Generate();
  } while(weight < random_generator.Rndm());
  std::vector<TLorentzVector> event(4);
  for(Int_t i = 0; i < 4; i++) {
    TLorentzVector *v = m_phasespace.GetDecay(i);
    event[i] = *v;
  }
  return event;
}
