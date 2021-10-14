// Martin Duy Tat 2nd November 2020

#include<vector>
#include<complex>
#include"Event.h"
#include"TLorentzVector.h"
#include"TMath.h"

Event::Event() {
  m_momenta = std::vector<double>(16, 0);
}

Event::Event(const std::vector<double> &p): m_momenta(p) {
}

Event::Event(const std::vector<double> &p, const std::complex<double> &D_amplitude, const std::complex<double> &DBAR_amplitude, double weight): m_momenta(p), m_Damplitude(D_amplitude), m_DBARamplitude(DBAR_amplitude), m_weight(weight) {
}

const std::vector<double>& Event::GetEventVector() const {
  return m_momenta;
}

Event::Event(const std::vector<TLorentzVector> &p) {
  m_momenta = std::vector<double>(16);
  for(int i = 0; i < 4; i++) {
    m_momenta[4*i + 0] = p[i][0];
    m_momenta[4*i + 1] = p[i][1];
    m_momenta[4*i + 2] = p[i][2];
    m_momenta[4*i + 3] = p[i][3];
  }
}

void Event::GetAmplitudes(std::complex<double> &D_amplitude, std::complex<double> &DBAR_amplitude) const {
  D_amplitude = m_Damplitude;
  DBAR_amplitude = m_DBARamplitude;
}

double Event::GetInvMass2(const int &particle1, const int &particle2) const {
  return TMath::Sqrt(TMath::Power(m_momenta[4*particle1 + 3] + m_momenta[4*particle2 + 3], 2)
		   - TMath::Power(m_momenta[4*particle1 + 0] + m_momenta[4*particle2 + 0], 2)
		   - TMath::Power(m_momenta[4*particle1 + 1] + m_momenta[4*particle2 + 1], 2)
		   - TMath::Power(m_momenta[4*particle1 + 2] + m_momenta[4*particle2 + 2], 2));
}

double Event::GetInvMass3(const int &particle1, const int &particle2, const int &particle3) const {
  return TMath::Sqrt(TMath::Power(m_momenta[4*particle1 + 3] + m_momenta[4*particle2 + 3] + m_momenta[4*particle3 + 3], 2)
		   - TMath::Power(m_momenta[4*particle1 + 0] + m_momenta[4*particle2 + 0] + m_momenta[4*particle3 + 0], 2)
		   - TMath::Power(m_momenta[4*particle1 + 1] + m_momenta[4*particle2 + 1] + m_momenta[4*particle3 + 1], 2)
		   - TMath::Power(m_momenta[4*particle1 + 2] + m_momenta[4*particle2 + 2] + m_momenta[4*particle3 + 2], 2));
}

double Event::GetWeight() const {
  return m_weight;
}
