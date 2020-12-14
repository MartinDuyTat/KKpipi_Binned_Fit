// Martin Duy Tat 2nd November 2020

#include"CPParameters.h"
#include"TMatrixD.h"

CPParameters::CPParameters(const double &xplus, const double &xminus, const double &yplus, const double &yminus): m_xplus(xplus), m_xminus(xminus), m_yplus(yplus), m_yminus(yminus) {
}

void CPParameters::GetCPParameters(double &xplus, double &xminus, double &yplus, double &yminus) const {
  xplus = m_xplus;
  xminus = m_xminus;
  yplus = m_yplus;
  yminus = m_yminus;
}

void CPParameters::SetCov(double *CovMatrix) {
  m_CovMatrix.ResizeTo(4, 4);
  m_CovMatrix.SetMatrixArray(CovMatrix);
}

TMatrixD CPParameters::GetCov() const {
  return m_CovMatrix;
}
