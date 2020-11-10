// Martin Duy Tat 10th November 2020

#include"Gamma.h"
#include"CPParameters.h"
#include"TMatrixD.h"

Gamma::Gamma(double rB, double deltaB, double gamma): m_rB(rB), m_deltaB(deltaB), m_gamma(gamma) {
}

void Gamma::GetGammaParameters(double &rB, double &deltaB, double &gamma) const {
  rB = m_rB;
  deltaB = m_deltaB;
  gamma = m_gamma;
}

void Gamma::SetCov(double *CovMatrix) {
  m_CovMatrix.ResizeTo(3, 3);
  m_CovMatrix.SetMatrixArray(CovMatrix);
}

TMatrixD Gamma::GetCov() const {
  return m_CovMatrix;
}
