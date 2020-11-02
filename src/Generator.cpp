// Martin Duy Tat

#include<vector>
#include"Generator.h"
#include"TGenPhaseSpace.h"
#include"TLorentzVector.h"

Generator::Generator(const Double_t &mass_parent, const Double_t *mass_decay, Int_t particles) {
  m_phasespace = TGenPhaseSpace();
  TLorentzVector P(mass_decay);
  m_phasespace.SetDecay(P, particles, mass_parent);
}

std::vector<TLorentzVector> Generator::Generator() {
  Double_t = weight;
  TRandom3 random_generator(0);
  do {
    weight = m_phasespace.Generate();
  } while(weight < random_generator.Rndm());
  std::vector<TLorentzVector> event(4);
  for(Int_t i = 0; i < 4; i++) {
    event[i] = m_phasespace.GetDecay(i);
  }
  return event;
}
