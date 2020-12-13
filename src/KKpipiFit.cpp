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
      phasespace = new AmplitudePhaseSpace(bins);
    } else {
      phasespace = nullptr;
    }
    return phasespace;
  }

  void LoadInputDataIntoBins(const std::string &Bplusfile, const std::string &Bminusfile, BinList &binlist) {
    TFile fBplus(Bplusfile.c_str(), "READ");
    TFile fBminus(Bminusfile.c_str(), "READ");
    TTree *treeBplus, *treeBminus;
    fBplus.GetObject("DalitzEventList", treeBplus);
    fBminus.GetObject("DalitzEventList", treeBminus);
    binlist.LoadTTree(treeBplus, +1);
    binlist.LoadTTree(treeBminus, -1);
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

  void SplitTree(TTree *tree, TTree *treeSmall, const int &StartEvent, const int &SampleSize) {
    std::vector<Double_t> four_momentum(16);
    Double_t *p = four_momentum.data();
    for(int i = 0; i < 16; i++) {
      std::string address = tree->GetListOfBranches()[0][i]->GetName();
      if(i%4 != 0) {
	tree->SetBranchAddress(address.c_str(), p + i - 1);
	treeSmall->Branch(address.c_str(), p + i - 1, (address + "/D").c_str());
      } else {
	tree->SetBranchAddress(address.c_str(), p + i + 3);
	treeSmall->Branch(address.c_str(), p + i + 3, (address + "/D").c_str());
      }
    }
    for(int i = StartEvent; i < StartEvent + SampleSize; i++) {
      tree->GetEntry(i);
      treeSmall->Fill();
    }
  }

}
