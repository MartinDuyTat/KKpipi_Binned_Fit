// Martin Duy Tat 2nd November 2020

#include"CPParameters.h"

CPParameters::CPParameters(double xplus, double xminus, double yplus, double yminus): m_xplus(xplus), m_xminus(xminus), m_yplus(yplus), m_yminus(yminus) {
}

void CPParameters::GetCPParameters(double &xplus, double &xminus, double &yplus, double &yminus) const {
  xplus = m_xplus;
  xminus = m_xminus;
  yplus = m_yplus;
  yminus = m_yminus;
}
