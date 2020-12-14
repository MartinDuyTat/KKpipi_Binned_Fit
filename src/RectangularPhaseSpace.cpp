// Martin Duy Tat 12th November 2020

#include<vector>
#include<numeric>
#include<algorithm>
#include<map>
#include"RectangularPhaseSpace.h"
#include"PhaseSpaceParameterisation.h"
#include"TLorentzVector.h"
#include"Event.h"
#include"TMath.h"
#include"TVector3.h"
#include"KKpipiMath.h"
#include"Constants.h"

RectangularPhaseSpace::RectangularPhaseSpace(const std::vector<int> &bins): PhaseSpaceParameterisation(std::accumulate(bins.begin(), bins.end(), 1, std::multiplies<int>())), m_Binning(bins) {
  m_xlow = {KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI, KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI, -1.0, -1.0, -TMath::Pi()};
  m_xhigh = {KKpipi_Constants::MASS_D - KKpipi_Constants::MASS_K - KKpipi_Constants::MASS_PI, KKpipi_Constants::MASS_D - KKpipi_Constants::MASS_K - KKpipi_Constants::MASS_PI, 1.0, 1.0, TMath::Pi()};
  m_BinMap.resize(5);
  std::vector<double> BinWidths(5);
  for(int i = 0; i < 5; i++) {
    BinWidths[i] = (m_xhigh[i] - m_xlow[i])/m_Binning[i];
  }
  for(int i = 0; i < 5; i++) {
    for(int j = 0; j < m_Binning[i]; j++) {
      m_BinMap[i].insert({m_xlow[i] + BinWidths[i]*(j + 1), j});
    }
  }
}

RectangularPhaseSpace::RectangularPhaseSpace(): RectangularPhaseSpace({1, 1, 1, 1, 1}) {
}

RectangularPhaseSpace::~RectangularPhaseSpace() {
}
  
int RectangularPhaseSpace::WhichBin(const Event &event) const {
  std::vector<double> X = KKpipiMath::RectCoordinates(event);
  return m_BinMap[0].lower_bound(X[0])->second + m_BinMap[1].lower_bound(X[1])->second*m_Binning[0] + m_BinMap[2].lower_bound(X[2])->second*m_Binning[0]*m_Binning[1] + m_BinMap[3].lower_bound(X[3])->second*m_Binning[0]*m_Binning[1]*m_Binning[2] + m_BinMap[4].lower_bound(X[4])->second*m_Binning[0]*m_Binning[1]*m_Binning[2]*m_Binning[3];
}

int RectangularPhaseSpace::NumberOfBins() const {
  return PhaseSpaceParameterisation::NumberOfBins();
}

double RectangularPhaseSpace::GetLowerBoundary(int coordinate) const {
  return m_xlow[coordinate];
}

double RectangularPhaseSpace::GetUpperBoundary(int coordinate) const {
  return m_xhigh[coordinate];
}
