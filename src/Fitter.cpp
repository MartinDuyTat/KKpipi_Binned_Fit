// Martin Duy Tat 2nd November 2020

#include"Fitter.h"
#include"BinList.h"
#include"DDecayParameters.h"
#include"CPParameters.h"
#include"Likelihood.h"
#include"Minuit2/Minuit2Minimizer.h"
#include"Math/Functor.h"

Fitter::Fitter(const BinList &binlist, const DDecayParameters &ddparameters): m_binlist(binlist), m_ddparameters(ddparameters) {
}

void Fitter::DoFit(CPParameters &cpparameters) {
  ROOT::Minuit2::Minuit2Minimizer *mini = new ROOT::Minuit2::Minuit2Minimizer();
  Likelihood *likelihood = new Likelihood(m_binlist, m_ddparameters);
  ROOT::Math::Functor fcn(*likelihood, 4);
  mini->SetFunction(fcn);
  double xplus, xminus, yplus, yminus;
  cpparameters.GetCPParameters(xplus, xminus, yplus, yminus);
  mini->SetVariable(0, "xplus", xplus, 1);
  mini->SetVariable(1, "xminus", xminus, 1);
  mini->SetVariable(2, "yplus", yplus, 1);
  mini->SetVariable(3, "yminus", yminus, 1);
  mini->Minimize();
  const double *xs = mini->X();
  cpparameters = CPParameters(xs[0], xs[1], xs[2], xs[3]);
  double CovMatrix[16];
  mini->GetCovMatrix(CovMatrix);
  cpparameters.SetCov(CovMatrix);
  delete mini;
  delete likelihood;
}
