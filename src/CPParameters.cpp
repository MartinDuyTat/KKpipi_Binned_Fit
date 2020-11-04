// Martin Duy Tat 2nd November 2020

#include"CPParameters.h"

CPParameters::CPParameters(double xplus, double xminus, double yplus, double yminus): m_xplus(xplus), m_xminus(xminus), m_yplus(yplus), m_yminus(yminus), m_xpluserror(0), m_xminuserror(0), m_ypluserror(0), m_yminuserror(0) {
}

void CPParameters::GetCPParameters(double &xplus, double &xminus, double &yplus, double &yminus) const {
  xplus = m_xplus;
  xminus = m_xminus;
  yplus = m_yplus;
  yminus = m_yminus;
}

void CPParameters::SetError(double xplus, double xminus, double yplus, double yminus) {
  m_xpluserror = xplus;
  m_xminuserror = xminus;
  m_ypluserror = yplus;
  m_yminuserror = yminus;
}

void CPParameters::GetError(double &xplus, double &xminus, double &yplus, double &yminus) const {
  xplus = m_xpluserror;
  xminus = m_xminuserror;
  yplus = m_ypluserror;
  yminus = m_yminuserror;
}
