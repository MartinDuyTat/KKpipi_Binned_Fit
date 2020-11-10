// Martin Duy Tat 10th November 2020

#include"FitGamma.h"
#include"Gamma.h"
#include"XYLikelihood.h"
#include"CPParameters.h"

FitGamma::FitGamma(CPParameters cpparameters): m_cpparameters(cpparameters) {
}

void FitGamma::DoFit(Gamma &GammaParams) {
  if(m_minimiser == nullptr) {
    m_minimiser = new ROOT::Minuit2::Minuit2Minimizer();
  }
  if(m_likelihood == nullptr) {
    m_likelihood = new XYLikelihood(m_cpparameters);
  }
  m_fcn = ROOT::Math::Functor(*m_likelihood, 3);
  m_minimiser->SetFunction(m_fcn);
  double rB, deltaB, gamma;
  GammaParams.GetGammaParameters(rB, deltaB, gamma);
  m_minimiser->SetVariable(0, "r_B", rB, 1);
  m_minimiser->SetVariable(1, "delta_B", deltaB, 10);
  m_minimiser->SetVariable(2, "gamma", gamma, 10);
  m_minimiser->Minimize();
  const double *result = m_minimiser->X();
  GammaParams = Gamma(result[0], result[1], result[2]);
  double CovMatrix[16];
  m_minimiser->GetCovMatrix(CovMatrix);
  GammaParams.SetCov(CovMatrix);
}

FitGamma::~FitGamma() {
  delete m_likelihood;
  m_likelihood = nullptr;
  delete m_minimiser;
  m_minimiser = nullptr;
}
