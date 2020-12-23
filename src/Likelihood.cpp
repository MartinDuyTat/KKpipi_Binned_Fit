// Martin Duy Tat 2nd November 2020

#include"Likelihood.h"
#include"DDecayParameters.h"
#include"CPParameters.h"
#include"BinList.h"
#include"TMath.h"
#include"Math/PdfFuncMathCore.h"
#include"KKpipiMath.h"

Likelihood::Likelihood(const BinList &bins, const DDecayParameters &ddparameters): m_bins(bins), m_ddparameters(ddparameters), m_leastsquares(false) {
}

double Likelihood::LogPoissonPDF(const int &x, const double &mu) const {
  return x*TMath::Log(mu) - ROOT::Math::lgamma(x + 1) - mu;
}

double Likelihood::operator()(const double *cpparameters) {
  double loglikelihood = 0;
  std::vector<int> eventsBplus = m_bins.GetEvents(+1, +1);
  std::vector<int> eventsBminus = m_bins.GetEvents(-1, +1);
  std::vector<int> eventsCPBplus = m_bins.GetEvents(+1, -1);
  std::vector<int> eventsCPBminus = m_bins.GetEvents(-1, -1);
  std::vector<double> predictedBplus, predictedBminus, predictedCPBplus, predictedCPBminus;
  int totalBplus = 0, totalBminus = 0;
  for(int i = 0; i < m_bins.NumberBins(); i++) {
    totalBplus += eventsBplus[i] + eventsCPBplus[i];
    totalBminus += eventsBminus[i] + eventsCPBminus[i];
  }
  CPParameters cpparam(cpparameters[0], cpparameters[1], cpparameters[2], cpparameters[3]);
  KKpipiMath::ExpectedNumberOfEvents(m_ddparameters, cpparam, totalBplus, totalBminus, predictedBplus, predictedBminus, predictedCPBplus, predictedCPBminus);
  if(m_leastsquares) {
    for(int i = 0; i < m_bins.NumberBins(); i++) {
      loglikelihood += TMath::Power(eventsBplus[i] - predictedBplus[i], 2)/predictedBplus[i];
      loglikelihood += TMath::Power(eventsBminus[i] - predictedBminus[i], 2)/predictedBminus[i];
      loglikelihood += TMath::Power(eventsCPBplus[i] - predictedCPBplus[i], 2)/predictedCPBplus[i];
      loglikelihood += TMath::Power(eventsCPBminus[i] - predictedCPBminus[i], 2)/predictedCPBminus[i];
    }
  } else {
    for(int i = 0; i < m_bins.NumberBins(); i++) {
      loglikelihood += -2*LogPoissonPDF(eventsBplus[i], predictedBplus[i]);
      loglikelihood += -2*LogPoissonPDF(eventsBminus[i], predictedBminus[i]);
      loglikelihood += -2*LogPoissonPDF(eventsCPBplus[i], predictedCPBplus[i]);
      loglikelihood += -2*LogPoissonPDF(eventsCPBminus[i], predictedCPBminus[i]);
    }
  }
  return loglikelihood;
}

void Likelihood::SetLeastSquares(bool UseLeastSquares) {
  m_leastsquares = UseLeastSquares;
}
