// Martin Duy Tat 25th November 2020

#include<vector>
#include<string>
#include<complex>
#include<fstream>
#include<stdlib.h>
#include"SophisticatedPhaseSpace.h"
#include"Event.h"
#include"Amplitude.h"
#include"TFile.h"
#include"TTree.h"
#include"KKpipiMath.h"
#include"TMath.h"

SophisticatedPhaseSpace::SophisticatedPhaseSpace(const int &nbins): PhaseSpaceParameterisation(nbins), m_regions(95), m_binregion(2) {
}

void SophisticatedPhaseSpace::ReadAverageStrongPhases(const std::string &filename) {
  std::ifstream f(filename);
  std::string line;
  std::getline(f, line);
  int N = atoi(line.c_str());
  m_LookupBins = std::vector<std::vector<std::vector<std::vector<int>>>>(m_regions, std::vector<std::vector<std::vector<int>>>(N, std::vector<std::vector<int>>(N, std::vector<int>(N))));
  std::vector<std::string> lines(m_regions + 3);
  double BinWidth = 2*TMath::Pi()/(NumberOfBins() - m_binregion);
  /*std::map<double, int> BinMap;
  for(int i = 0; i < NumberOfBins() - m_binregion; i++) {
    BinMap.insert({-TMath::Pi() + BinWidth*(i + 1), i});
    }*/
  bool endfile = true;
  while(endfile) {
    for(int i = 0; i < m_regions + 2; i++ ) {
      endfile = endfile && std::getline(f, lines[i], ',');
    }
    endfile = endfile && std::getline(f, lines[m_regions + 2]);
    for(int i = 0; i < m_regions; i++) {
      double phase = atof(lines[i + 3].c_str());
      //m_LookupBins[i][atoi(lines[0].c_str())][atoi(lines[1].c_str())][atoi(lines[2].c_str())] = BinMap.lower_bound(phase)->second;
      m_LookupBins[i][atoi(lines[0].c_str())][atoi(lines[1].c_str())][atoi(lines[2].c_str())] = static_cast<int>((phase + TMath::Pi())/BinWidth);
    }
  }
}

SophisticatedPhaseSpace::SophisticatedPhaseSpace(const int &nbins, const std::string &filename): SophisticatedPhaseSpace(nbins) {
  ReadAverageStrongPhases(filename);
}

SophisticatedPhaseSpace::~SophisticatedPhaseSpace() {
  if(m_AmplitudeModel != nullptr) {
    delete m_AmplitudeModel;
  }
}

void SophisticatedPhaseSpace::x3x4WhichBinGrid(const Event &event, int &x3bin, int &x4bin) const {
  std::vector<double> X = KKpipiMath::RectCoordinates(event.GetEventVector());
  double x3lower = -1.0;
  double x3upper = 1.0;
  double x4lower = -1.0;
  double x4upper = 1.0;
  x3bin = int((X[2] - x3lower)*m_x3bins/(x3upper - x3lower));
  x4bin = int((X[3] - x4lower)*m_x4bins/(x4upper - x4lower));
}

void SophisticatedPhaseSpace::LoadAmplitudeModel() {
  m_AmplitudeModel = new Amplitude;
}

double SophisticatedPhaseSpace::StrongPhase(const Event &event) const {
  std::vector<double> fourmomenta = event.GetEventVector();
  return std::arg(m_AmplitudeModel->operator()(fourmomenta, +1)) - std::arg(m_AmplitudeModel->operator()(fourmomenta, -1));
}

void SophisticatedPhaseSpace::CalculateStrongPhases(const std::string &BplusFilename, const std::string &BminusFilename, const std::string &MeanFilename, const std::string &RMSFilename) const {
  std::vector<std::vector<double>> average(m_x3bins, std::vector<double>(m_x4bins));
  std::vector<std::vector<double>> rms(m_x3bins, std::vector<double>(m_x4bins));
  std::vector<std::vector<int>> numberevents(m_x3bins, std::vector<int>(m_x4bins));
  TFile fBplus(BplusFilename.c_str(), "READ");
  TFile fBminus(BminusFilename.c_str(), "READ");
  TTree *treeBplus, *treeBminus;
  fBplus.GetObject("DalitzEventList", treeBplus);
  fBminus.GetObject("DalitzEventList", treeBminus);
  std::vector<double> FourMomentumBplus(16), FourMomentumBminus(16);
  for(int i = 0; i < 16; i++) {
    std::string address = treeBplus->GetListOfBranches()[0][i]->GetName();
    if(i%4 == 0) {
      treeBplus->SetBranchAddress(address.c_str(), FourMomentumBplus.data() + i + 3);
      treeBminus->SetBranchAddress(address.c_str(), FourMomentumBminus.data() + i + 3);
    } else {
      treeBplus->SetBranchAddress(address.c_str(), FourMomentumBplus.data() + i - 1);
      treeBminus->SetBranchAddress(address.c_str(), FourMomentumBminus.data() + i - 1);
    }
  }
  int x3bin, x4bin;
  double BplusStrongPhase, BminusStrongPhase;
  for(int i = 0; i < treeBplus->GetEntries(); i++) {
    treeBplus->GetEvent(i);
    treeBminus->GetEvent(i);
    Event EventBplus(FourMomentumBplus), EventBminus(FourMomentumBminus);
    BplusStrongPhase = StrongPhase(EventBplus);
    BminusStrongPhase = StrongPhase(EventBminus);
    x3x4WhichBinGrid(EventBplus, x3bin, x4bin);
    numberevents[x3bin][x4bin] += 1;
    average[x3bin][x4bin] += BplusStrongPhase;
    rms[x3bin][x4bin] += BplusStrongPhase*BplusStrongPhase;
    x3x4WhichBinGrid(EventBminus, x3bin, x4bin);
    numberevents[x3bin][x4bin] += 1;
    average[x3bin][x4bin] += BminusStrongPhase;
    rms[x3bin][x4bin] += BminusStrongPhase*BminusStrongPhase;
  }
  std::ofstream AverageFile(MeanFilename), RMSFile(RMSFilename);
  for(int i = 0; i < m_x3bins; i++) {
    for(int j = 0; j < m_x4bins; j++) {
      average[i][j] = average[i][j]/numberevents[i][j];
      rms[i][j] = (rms[i][j]/numberevents[i][j]) - average[i][j]*average[i][j];
      AverageFile << average[i][j] << " ";
      RMSFile << rms[i][j] << " ";
    }
    AverageFile << std::endl;
    RMSFile << std::endl;
  }
  AverageFile.close();
  RMSFile.close();
}

void SophisticatedPhaseSpace::ClearAverageStrongPhases() {
  std::vector<std::vector<std::vector<std::vector<int>>>>().swap(m_LookupBins);
}

int SophisticatedPhaseSpace::NumberOfRegions() const {
  return m_regions;
}

/*int SophisticatedPhaseSpace::WhichRegion(const std::vector<double> &X) const {
  if(X[2] + X[3] > 1.4) {
    return 0;
  } else if(10*X[2] - 7*X[3] < -10) {
    return 1;
  } else if(45*X[2] - 35*X[3] < -17) { 
    return 2;
  } else if(7*X[2] - 5*X[3] > 5) {
    return 3;
  } else if(9*X[2] + 4*X[3] > 5) {
    return 4;
  } else {
    return 5;
  }
}*/

int SophisticatedPhaseSpace::WhichRegion(const std::vector<double> &X) const {
  if(X[2] > 0.4 && X[3] > 0.4) {
    if(X[2] + X[3] > 1.4) {
      if(X[2] > X[3]) {
	return 0;
      } else {
	return 1;
      }
    } else {
      if(X[2] > X[3]) {
	return 2;
      } else {
	return 3;
      }
    }
  }
  int x3bin, x4bin;
  if(X[3] < 0.4) {
    x3bin = static_cast<int>((X[2] + 1.0)/(2.0/10.0));
    x4bin = static_cast<int>((X[3] + 1.0)/(1.4/7.0));
    return x3bin + 10*x4bin + 4;
  } else {
    x3bin = static_cast<int>((X[2] + 1.0)/(1.4/7.0));
    x4bin = static_cast<int>((X[3] - 0.4)/(0.6/3.0));
    return x3bin + 7*x4bin + 74;
  }
}

int SophisticatedPhaseSpace::WhichBin(const Event &event) const {
  std::vector<double> X = KKpipiMath::RectCoordinates(event);
  int Region = WhichRegion(X);
  if(Region < m_binregion) {
    return Region;
  }
  int N = m_LookupBins[0].size();
  double dx1 = (RectangularPhaseSpace::GetUpperBoundary(0) - RectangularPhaseSpace::GetLowerBoundary(0))/N;
  double dx2 = (RectangularPhaseSpace::GetUpperBoundary(1) - RectangularPhaseSpace::GetLowerBoundary(1))/N;
  double dx5 = (RectangularPhaseSpace::GetUpperBoundary(4) - RectangularPhaseSpace::GetLowerBoundary(4))/N;
  int x1_bin = static_cast<int>((X[0] - RectangularPhaseSpace::GetLowerBoundary(0))/dx1);
  int x2_bin = static_cast<int>((X[1] - RectangularPhaseSpace::GetLowerBoundary(1))/dx2);
  int x5_bin = static_cast<int>((X[4] - RectangularPhaseSpace::GetLowerBoundary(4))/dx5);
  return m_LookupBins[Region][x1_bin][x2_bin][x5_bin] + m_binregion;
}
  
  

int SophisticatedPhaseSpace::NumberOfBins() const {
  return PhaseSpaceParameterisation::NumberOfBins();
}
