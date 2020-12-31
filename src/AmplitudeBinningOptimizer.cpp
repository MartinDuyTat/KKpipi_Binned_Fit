// Martin Duy Tat 31st December 2020

#include<string>
#include<vector>
#include"TMath.h"
#include"Minuit2/Minuit2Minimizer.h"
#include"Math/Functor.h"
#include"KKpipiMath.h"
#include"AmplitudeBinningOptimizer.h"
#include"AmplitudePhaseSpace.h"
#include"PhaseSpaceParameterisation.h"
#include"DDecayParameters.h"

AmplitudeBinningOptimizer::AmplitudeBinningOptimizer(const int &bins, const std::string filename, const std::vector<double> &BinEdges): m_bins(bins) {
  double BinWidth = 2*TMath::Pi()/m_bins;
  if(!BinEdges.empty()) {
    m_BinEdges = BinEdges;
  } else {
    for(int i = 0; i < (m_bins - 2)/2; i++) {
      m_BinEdges.push_back((i + 1)*BinWidth);
    }
  }
  m_eventlist.LoadTree(filename);
}

double AmplitudeBinningOptimizer::operator ()(const double *BinEdges) {
  AmplitudePhaseSpace aph(m_bins);
  aph.ReadAmplitudeFromEvent(true);
  aph.SetBinEdges(std::vector<double>(BinEdges, BinEdges + (m_bins - 2)/2));
  aph.UseVariableBinWidths(true);
  PhaseSpaceParameterisation *psp = &aph;
  DDecayParameters ddparameters(psp, m_eventlist);
  return -KKpipiMath::CalculateExactBinningQValue(ddparameters);
}

double AmplitudeBinningOptimizer::OptimizeBinEdges() {
  ROOT::Minuit2::Minuit2Minimizer *mini = new ROOT::Minuit2::Minuit2Minimizer();
  ROOT::Math::Functor fcn(*this, (m_bins - 2)/2);
  mini->SetFunction(fcn);
  for(int i = 0; i < (m_bins - 2)/2; i++) {
    mini->SetVariable(i, "EdgeParameter" + std::to_string(i), m_BinEdges[i], 2*TMath::Pi()/m_bins);
  }
  mini->Minimize();
  const double *xs = mini->X();
  m_BinEdges = std::vector<double>(xs, xs + (m_bins - 2)/2);
  double MinValue = mini->MinValue();
  delete mini;
  return -MinValue;
}

const std::vector<double>& AmplitudeBinningOptimizer::GetBinEdges() const {
  return m_BinEdges;
}
