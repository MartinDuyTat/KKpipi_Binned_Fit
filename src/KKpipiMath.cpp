// Martin Duy Tat 1st December 2020

#include<vector>
#include"KKpipiMath.h"
#include"TLorentzVector.h"
#include"Event.h"
#include"Constants.h"
#include"TVector3.h"
#include"TMath.h"

namespace KKpipiMath {

  std::vector<TLorentzVector> ConvertEventTo4Vectors(const std::vector<double> &event) {
    std::vector<TLorentzVector> momenta(4);
    for(int i = 0; i < 4; i++) {
      momenta[i] = TLorentzVector(event.data() + 4*i);
    }
    return momenta;
  }

  std::vector<double> RectCoordinates(const Event &event) {
    std::vector<double> X(5);
    std::vector<TLorentzVector> P = ConvertEventTo4Vectors(event.GetEventVector());
    // Use invariant mass of K+pi+ and K-pi- as variables
    double mplus = event.GetInvMass2(0, 2);
    double mminus = event.GetInvMass2(1, 3);
    double mmin = std::min(mplus, mminus) - KKpipi_Constants::MASS_K - KKpipi_Constants::MASS_PI;
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

  std::vector<double> ConvertXToMomenta(const std::vector<double> &X) {
    // Find the invariant masses of K+pi+ and K-pi-
    double alpha = (std::min(X[0], X[1]) - KKpipi_Constants::MASS_K - KKpipi_Constants::MASS_PI)/2.0;
    double mplus = X[0] - alpha;
    double mminus = X[1] - alpha;
    // Speed of D meson in the K+pi+ rest frame
    double betap = TMath::Sqrt(mplus*mplus - (KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI)*(KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI))/mplus;
    // Energy of K+ in the K+pi+ rest frame
    double EKp = (mplus*mplus + KKpipi_Constants::MASS_K*KKpipi_Constants::MASS_K - KKpipi_Constants::MASS_PI*KKpipi_Constants::MASS_PI)/(2*mplus);
    // Momentum of K+ in the K+pi+ rest frame
    double pKp = TMath::Sqrt(EKp*EKp - KKpipi_Constants::MASS_K*KKpipi_Constants::MASS_K);
    // Energy of pi+ in the K+pi+ rest frame
    double Epip = (mplus*mplus - KKpipi_Constants::MASS_K*KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI*KKpipi_Constants::MASS_PI)/(2*mplus);
    // Put into four-vectors, with K+ at an angle theta+ relative to the D meson, and choose them in the xz plane
    TLorentzVector P_Kp(pKp*TMath::Sqrt(1 - X[2]*X[2]), 0.0, -pKp*X[2], EKp);
    TLorentzVector P_pip(-pKp*TMath::Sqrt(1 - X[2]*X[2]), 0.0, pKp*X[2], Epip);
    // Speed of D meson in the K-pi- rest frame
    double betam = TMath::Sqrt(mminus*mminus - (KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI)*(KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI))/mminus;
    // Energy of K- in the K-pi- rest frame
    double EKm = (mminus*mminus + KKpipi_Constants::MASS_K*KKpipi_Constants::MASS_K - KKpipi_Constants::MASS_PI*KKpipi_Constants::MASS_PI)/(2*mminus);
    // Momentum of K- in the K-pi- rest frame
    double pKm = TMath::Sqrt(EKm*EKm - KKpipi_Constants::MASS_K*KKpipi_Constants::MASS_K);
    // Energy of pi- in the K-pi- rest frame
    double Epim = (mminus*mminus - KKpipi_Constants::MASS_K*KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI*KKpipi_Constants::MASS_PI)/(2*mminus);
    // Put into four-vector, with K- at an angle theta- relative to the D meson, and put them initially in the xz plane
    TLorentzVector P_Km(-pKm*TMath::Sqrt(1 - X[3]*X[3]), 0.0, pKm*X[3], EKm);
    TLorentzVector P_pim(pKm*TMath::Sqrt(1 - X[3]*X[3]), 0.0, -pKm*X[3], Epim);
    // Boost K+ and pi+ to the D meson rest frame
    P_Kp.Boost(0.0, 0.0, betap);
    P_pip.Boost(0.0, 0.0, betap);
    // Rotate the K- and pi- by phi because their decay planes are not parallel
    P_Km.RotateZ(-X[4]);
    P_pim.RotateZ(-X[4]);
    // Boost I- and pi- to the D meson rest frame
    P_Km.Boost(0.0, 0.0, -betam);
    P_pim.Boost(0.0, 0.0, -betam);
    // Put everything together and return it
    std::vector<double> event = {P_Kp[0], P_Kp[1], P_Kp[2], P_Kp[3], P_Km[0], P_Km[1], P_Km[2], P_Km[3], P_pip[0], P_pip[1], P_pip[2], P_pip[3], P_pim[0], P_pim[1], P_pim[2], P_pim[3]};
    return event;
  }
}
