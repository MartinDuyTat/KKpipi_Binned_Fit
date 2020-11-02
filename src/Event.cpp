// Martin Duy Tat 2nd November 2020

#include"Event.h"
#include"TLorentzVector.h"
#include"TMath.h"

Event::Event() {
  m_momenta = std::vector<double>(16, 0);
}

Event::Event(std::vector<double> p): m_momenta(p) {
}

std::vector<double> Event::GetEventVector() {
  return m_momenta;
}

Event::Event(const Event &event): m_momenta(event.GetEventVector()) {
}

Event::Event(const std::vector<TLorentzVector> &p) {
  for(int i = 0; i < 4; i++) {
    m_momenta[4*i + 0] = p[i][3];
    m_momenta[4*i + 1] = p[i][0];
    m_momenta[4*i + 2] = p[i][1];
    m_momenta[4*i + 3] = p[i][2];
  }
}

double Event::GetInvMass(int particle1, int particle2) {
  return TMath::Sqrt(TMath::Power(m_momenta[4*particle1 + 0] + m_momenta[4*particle2 + 0], 2) - TMath::Power(m_momenta[4*particle1 + 1] + m_momenta[4*particle2 + 1], 2) - TMath::Power(m_momenta[4*particle1 + 2] + m_momenta[4*particle2 + 2], 2) - TMath::Power(m_momenta[4*particle1 + 3] + m_momenta[4*particle2 + 3], 2));
}

double Event::GetInvMass(int particle1, int particle2) {
  return TMath::Sqrt(TMath::Power(m_momenta[4*particle1 + 0] + m_momenta[4*particle2 + 0] + m_momenta[4*particle3 + 0], 2) - TMath::Power(m_momenta[4*particle1 + 1] + m_momenta[4*particle2 + 1] + m_momenta[4*particle3 + 1], 2) - TMath::Power(m_momenta[4*particle1 + 2] + m_momenta[4*particle2 + 2] + m_momenta[4*particle3 + 2], 2) - TMath::Power(m_momenta[4*particle1 + 3] + m_momenta[4*particle2 + 3] + m_momenta[4*particle3 + 3], 2));
}
