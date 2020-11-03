// Martin Duy Tat 2nd November 2020

#include"Likelihood.h"
#include"DDecayParameters.h"
#include"CPParameters.h"
#include"BinList.h"
#include"TMath.h"
#include"Math/PdfFuncMathCore.h"

Likelihood::Likelihood(BinList bins, DDecayParameters ddparameters): m_bins(bins), m_ddparameters(ddparameters) {
}

double Likelihood::operator()(const double *cpparameters) {
  double loglikelihood = 0;
  std::vector<int> eventsBplus = m_bins.GetEvents(+1);
  std::vector<int> eventsBminus = m_bins.GetEvents(-1);
  std::vector<double> predictedBplus, predictedBminus;
  int totalBplus = 0, totalBminus = 0;
  for(int i = 0; i < m_bins.NumberBins(); i++) {
    totalBplus += eventsBplus[i];
    totalBminus += eventsBminus[i];
  }
  CPParameters cpparam(cpparameters[0], cpparameters[1], cpparameters[2], cpparameters[3]);
  m_bins.Predict(m_ddparameters, cpparam, predictedBplus, predictedBminus, totalBplus, totalBminus);
  for(int i = 0; i < m_bins.NumberBins(); i++) {
    loglikelihood += -2*TMath::Log(ROOT::Math::poisson_pdf(eventsBplus[i], predictedBplus[i]));
    loglikelihood += -2*TMath::Log(ROOT::Math::poisson_pdf(eventsBminus[i], predictedBminus[i]));
  }
  return loglikelihood;
}
