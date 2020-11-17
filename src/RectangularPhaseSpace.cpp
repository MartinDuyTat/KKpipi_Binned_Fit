// Martin Duy Tat 12th November 2020

#include<vector>
#include<numeric>
#include<algorithm>
#include<map>
#include"RectangularPhaseSpace.h"
#include"PhaseSpaceParameterisation.h"
#include"TLorentzVector.h"
#include"Event.h"
#include"TMath.h"
#include"TVector3.h"

RectangularPhaseSpace::RectangularPhaseSpace(std::vector<int> bins, double *masses): PhaseSpaceParameterisation(std::accumulate(bins.begin(), bins.end(), 1, std::multiplies<int>())), m_Binning(bins) {
  if(masses == nullptr) {
    m_Masses.resize(5);
    m_Masses = {1.86484, 0.493677, 0.493677, 0.13957018, 0.13957018};
  } else {
    m_Masses = std::vector<double>(masses, masses + 5);
  }
  m_xlow = {m_Masses[1] + m_Masses[3], m_Masses[2] + m_Masses[4], -1.0, -1.0, -TMath::Pi()};
  m_xhigh = {m_Masses[0] - m_Masses[2] - m_Masses[4], m_Masses[0] - m_Masses[1] - m_Masses[3], 1.0, 1.0, TMath::Pi()};
  m_BinMap.resize(5);
  std::vector<double> BinWidths(5);
  for(int i = 0; i < 5; i++) {
    BinWidths[i] = (m_xhigh[i] - m_xlow[i])/m_Binning[i];
  }
  for(int i = 0; i < 5; i++) {
    for(int j = 0; j < m_Binning[i]; j++) {
      m_BinMap[i].insert({m_xlow[i] + BinWidths[i]*(j + 1), j});
    }
  }
}

std::vector<TLorentzVector> RectangularPhaseSpace::ConvertTo4Vectors(const Event &event) const {
  std::vector<double> p = event.GetEvent();
  std::vector<TLorentzVector> momenta(4);
  for(int i = 0; i < 4; i++) {
    momenta[i] = TLorentzVector(p.data() + 4*i);
  }
  return momenta;
}

std::vector<double> RectangularPhaseSpace::RectCoordinates(const Event &event) const {
  std::vector<double> X(5);
  std::vector<TLorentzVector> P = ConvertTo4Vectors(event);
  // Use invariant mass of K+pi+ and K-pi- as variables
  double mplus = event.GetInvMass2(0, 2);
  double mminus = event.GetInvMass2(1, 3);
  double mmin = std::min(mplus, mminus) - m_Masses[1] - m_Masses[3];
  // Expand triangle in phase space into rectangle
  X[0] = mplus + mmin;
  X[1] = mminus + mmin;
  // 4-vector of D meson in the rest frame
  TLorentzVector Ptot = P[0] + P[1] + P[2] + P[3];
  // Get boost beta
  TVector3 beta = (P[0] + P[2]).BoostVector();
  TLorentzVector PtotTemp = Ptot, PKTemp = P[0];
  // Boost to rest frame of K+pi+ system
  PtotTemp.Boost(-beta);
  PKTemp.Boost(-beta);
  // Find angle between D meson and K+
  X[2] = PtotTemp.Vect().Unit().Dot(PKTemp.Vect().Unit());
  beta = (P[1] + P[3]).BoostVector();
  PtotTemp = Ptot;
  PKTemp = P[1];
  // Boost to rest frame of K-pi- system
  PtotTemp.Boost(-beta);
  PKTemp.Boost(-beta);
  // Find angle between D meson and K-
  X[3] = PtotTemp.Vect().Unit().Dot(PKTemp.Vect().Unit());
  // Find cosine and phi of the angle between decay planes of K+pi+ and K-pi-
  double sphi = P[0].Vect().Cross(P[2].Vect()).Unit().Cross(P[1].Vect().Cross(P[3].Vect()).Unit()).Dot((P[1] + P[3]).Vect().Unit());
  double cphi = P[0].Vect().Cross(P[2].Vect()).Unit().Dot(P[1].Vect().Cross(P[3].Vect()).Unit());
  X[4] = TMath::ATan2(sphi, cphi);
  return X;
}
  
int RectangularPhaseSpace::WhichBin(const Event &event) const {
  std::vector<double> X = RectCoordinates(event);
  return m_BinMap[0].lower_bound(X[0])->second + m_BinMap[1].lower_bound(X[1])->second*m_Binning[0] + m_BinMap[2].lower_bound(X[2])->second*m_Binning[0]*m_Binning[1] + m_BinMap[3].lower_bound(X[3])->second*m_Binning[0]*m_Binning[1]*m_Binning[2] + m_BinMap[4].lower_bound(X[4])->second*m_Binning[0]*m_Binning[1]*m_Binning[2]*m_Binning[3];
}

int RectangularPhaseSpace::NumberOfBins() const {
  return PhaseSpaceParameterisation::NumberOfBins();
}
