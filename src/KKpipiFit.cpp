// Martin Duy Tat 12th December 2020

#include<iostream>
#include<string>
#include<vector>
#include"TMath.h"
#include"TFile.h"
#include"TTree.h"
#include"PhaseSpaceParameterisation.h"
#include"NaivePhaseSpace.h"
#include"RectangularPhaseSpace.h"
#include"SophisticatedPhaseSpace.h"
#include"AmplitudePhaseSpace.h"
#include"BinList.h"
#include"CPParameters.h"
#include"Gamma.h"
#include"FitGamma.h"

namespace KKpipiFit {

  PhaseSpaceParameterisation* PickBinningScheme(const std::string &binning_choice) {
    PhaseSpaceParameterisation *phasespace;
    if(binning_choice == "naive") {
      std::cout << "Using Naive binning scheme\n";
      phasespace = new NaivePhaseSpace;
    } else if(binning_choice == "rectangular") {
      std::cout << "Using Rectangular binning scheme\n";
      std::vector<int> rectbins(5);
      std::cout << "Please input bins in each direction:\n";
      std::cin >> rectbins[0] >> rectbins[1] >> rectbins [2] >> rectbins[3] >> rectbins [4];
      phasespace = new RectangularPhaseSpace(rectbins);
    } else if(binning_choice == "sophisticated") {
      std::cout << "Using Sophisticated binning scheme\n";
      int bins;
      std::string phasefilename;
      std::cout << "Please input total number of bins:\n";
      std::cin >> bins;
      std::cout << "Please input file with mean strong phases:\n";
      std::cin >> phasefilename;
      phasespace = new SophisticatedPhaseSpace(bins, phasefilename);
    } else if(binning_choice == "amplitude") {
      std::cout << "Using Amplitude binning scheme\n";
      int bins;
      std::cout << "Please input total number of bins:\n";
      std::cin >> bins;
      std::string answer;
      std::cout << "Is binning along rD used?\n";
      std::cin >> answer;
      bool rDBinning = answer == "yes";
      phasespace = new AmplitudePhaseSpace(bins, rDBinning);
      AmplitudePhaseSpace *aph = static_cast<AmplitudePhaseSpace*>(phasespace);
      std::cout << "Read amplitudes from event?\n";
      std::cin >> answer;
      if(answer == "yes") {
	aph->ReadAmplitudeFromEvent(true);
      }
      std::cout << "Use variable bin widths?\n";
      std::cin >> answer;
      if(answer == "yes") {
	if(bins%2 != 0) {
	  std::cout << "Need an even number of bins for variable bin widths\n";
	  return nullptr;
	}
	std::vector<double> BinEdges((aph->NumberOfBins() - 2)/2);
	std::cout << "Please input the positive upper bin edges, in increasing order:\n";
	for(auto &edge : BinEdges) {
	  std::cin >> edge;
	}
	aph->SetBinEdges(BinEdges);
	aph->UseVariableBinWidths(true);
      }
    } else {
      phasespace = nullptr;
    }
    return phasespace;
  }

  void LoadTreesIntoBins(TTree *TreeBplus, TTree *TreeBminus, BinList &binlist, const int &StartEvent = 0, const int &TotalEvents = -1) {
    binlist.LoadTTree(TreeBplus, +1, StartEvent, TotalEvents);
    binlist.LoadTTree(TreeBminus, -1, StartEvent, TotalEvents);
    return;
  }

  void LoadInputDataIntoBins(const std::string &BplusFilename, const std::string &BminusFilename, BinList &binlist, const int &StartEvent, const int &TotalEvents) {
    TFile fBplus(BplusFilename.c_str(), "$READ");
    TFile fBminus(BminusFilename.c_str(), "READ");
    TTree *treeBplus = nullptr, *treeBminus = nullptr;
    fBplus.GetObject("DalitzEventList", treeBplus);
    fBminus.GetObject("DalitzEventList", treeBminus);
    LoadTreesIntoBins(treeBplus, treeBminus, binlist, StartEvent, TotalEvents);
    fBplus.Close();
    fBminus.Close();
    return;
  }

  void PrintXY(const CPParameters &cpparameters) {
    double xplus, xminus, yplus, yminus;
    cpparameters.GetCPParameters(xplus, xminus, yplus, yminus);
    TMatrixD cov = cpparameters.GetCov();
    std::cout << "Fitted x and y parameters:\n";
    std::cout << "xplus = " << xplus << " \u00B1 " << TMath::Sqrt(cov(0, 0)) << std::endl;
    std::cout << "xminus = " << xminus << " \u00B1 " << TMath::Sqrt(cov(1, 1)) << std::endl;
    std::cout << "yplus = " << yplus << " \u00B1 " << TMath::Sqrt(cov(2, 2)) << std::endl;
    std::cout << "yminus = " << yminus << " \u00B1 " << TMath::Sqrt(cov(3, 3)) << std::endl;
    return;
  }

  void PrintGamma (const Gamma &gammaparams) {
    double rB, deltaB, gamma;
    gammaparams.GetGammaParameters(rB, deltaB, gamma);
    TMatrixD gammacov = gammaparams.GetCov();
    std::cout << "Fitted r_B, delta_B and gamma parameters:\n";
    std::cout << "r_B = " << rB << " \u00B1 " << TMath::Sqrt(gammacov(0, 0)) << std::endl;
    std::cout << "delta_B = " << deltaB << " \u00B1 " << TMath::Sqrt(gammacov(1, 1)) << std::endl;
    std::cout << "gamma = " << gamma << " \u00B1 " << TMath::Sqrt(gammacov(2, 2)) << std::endl;
    return;
  }

  void DrawContours(const FitGamma &fitgamma) {
    std::string drawcontours;
    std::cout << "Draw contours?\n";
    std::cin >> drawcontours;
    if(drawcontours == "yes") {
      std::cout << "Drawing contours...\n";
      fitgamma.PlotContours("Contour_rB_vs_dB.png", "Contour_dB_vs_gamma.png", "Contour_gamma_vs_rB.png", 20);
      std::cout << "Finished drawing contours\n";
    }
    return;
  }

}
