// Martin Duy Tat 10th November 2020, based on code from Sneha Malde

#include"XYLikelihood.h"
#include"CPParameters.h"
#include"TMatrixD.h"
#include"TMath.h"

XYLikelihood::XYLikelihood(CPParameters cpparameters) {
  cpparameters.GetCPParameters(m_xplus, m_xminus, m_yplus, m_yminus);
  m_InvCov.ResizeTo(4, 4);
  m_InvCov = cpparameters.GetCov().Invert();
}

double XYLikelihood::operator()(const double *gamma) {
  TMatrixD diff(4, 1);
  diff(0, 0) = gamma[0]*TMath::Cos((gamma[1] + gamma[2])*TMath::Pi()/180.0) - m_xplus;
  diff(1, 0) = gamma[0]*TMath::Cos((gamma[1] - gamma[2])*TMath::Pi()/180.0) - m_xminus;
  diff(2, 0) = gamma[0]*TMath::Sin((gamma[1] + gamma[2])*TMath::Pi()/180.0) - m_yplus;
  diff(3, 0) = gamma[0]*TMath::Sin((gamma[1] - gamma[2])*TMath::Pi()/180.0) - m_yminus;
  TMatrixD diffT = diff;  
  TMatrixD Chi2 = diffT.T()*m_InvCov*diff;
  return Chi2.GetMatrixArray()[0];
}
