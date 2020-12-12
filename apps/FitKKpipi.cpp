// Martin Duy Tat 11th December 2020
/**
 * FitAmplitude is the program for doing the binned fitting using the amplitude model directly for binning
 * D meson decay parameters are obtained from an input file
 * @param 1 Binning scheme to use, in lowercase letters
 * @param 2 Filename of B+ event file
 * @param 3 Filename of B- event file
 * @param 4 Filename of D meson hadronic decay parameters
 */

#include<string>
#include<iostream>
#include<stdlib.h>
#include"KKpipiFit.h"
#include"PhaseSpaceParameterisation.h"
#include"NaivePhaseSpace.h"
#include"RectangularPhaseSpace.h"
#include"SophisticatedPhaseSpace.h"
#include"AmplitudePhaseSpace.h"
#include"TFile.h"
#include"TTree.h"
#include"TMath.h"
#include"TMatrixD.h"
#include"BinList.h"
#include"DDecayParameters.h"
#include"CPParameters.h"
#include"Fitter.h"
#include"FitGamma.h"
#include"Gamma.h"
#include"AmplitudePhaseSpace.h"

int main(int argc, char *argv[]) {
  if(argc != 5) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  std::cout << "Starting B->DK, D->KKpipi binned fit\n";
  std::cout << "Loading binning scheme for phase space...\n";
  PhaseSpaceParameterisation *phasespace = KKpipiFit::PickBinningScheme(std::string(argv[1]));
  if(phasespace == nullptr) {
    std::cout << "Invalid binning scheme\n";
    return 0;
  }
  std::cout << "Binning scheme loaded\n";
  std::cout << "Loading input data...\n";
  std::string Bplusfile = argv[2];
  std::string Bminusfile = argv[3];
  TFile fBplus(Bplusfile.c_str(), "READ");
  TFile fBminus(Bminusfile.c_str(), "READ");
  std::cout << "Input data ready\n";
  std::cout << "Putting input events into bins...\n";
  TTree *treeBplus, *treeBminus;
  fBplus.GetObject("DalitzEventList", treeBplus);
  fBminus.GetObject("DalitzEventList", treeBminus);
  BinList binlist(phasespace);
  binlist.LoadTTree(treeBplus, +1);
  binlist.LoadTTree(treeBminus, -1);
  std::cout << "Input events sorted into bins\n";
  std::cout << "Loading hadronic D decay parameters...\n";
  DDecayParameters ddparameters(argv[4]);
  std::cout << "Hadronic D parameters loaded\n";
  std::cout << "Loading fitter...\n";
  Fitter fit(binlist, ddparameters);
  CPParameters cpparameters(-0.09, 0.06, -0.04, 0.08);
  double xplus, xminus, yplus, yminus;
  std::cout << "Ready for fitting\n";
  std::cout << "Fitting x and y...\n";
  fit.DoFit(cpparameters);
  std::cout << "x and y fitted\n";
  cpparameters.GetCPParameters(xplus, xminus, yplus, yminus);
  TMatrixD cov = cpparameters.GetCov();
  std::cout << "Fitted x and y parameters:\n";
  std::cout << "xplus = " << xplus << " \u00B1 " << TMath::Sqrt(cov(0, 0)) << std::endl;
  std::cout << "xminus = " << xminus << " \u00B1 " << TMath::Sqrt(cov(1, 1)) << std::endl;
  std::cout << "yplus = " << yplus << " \u00B1 " << TMath::Sqrt(cov(2, 2)) << std::endl;
  std::cout << "yminus = " << yminus << " \u00B1 " << TMath::Sqrt(cov(3, 3)) << std::endl;
  std::cout << "Fitting r_B, delta_B and gamma...\n";
  FitGamma fitgamma(cpparameters);
  Gamma gammaparams(0.1, 130.0, 75.0);
  fitgamma.DoFit(gammaparams);
  std::cout << "r_B, delta_B and gamma fitted\n";
  double rB, deltaB, gamma;
  gammaparams.GetGammaParameters(rB, deltaB, gamma);
  TMatrixD gammacov = gammaparams.GetCov();
  std::cout << "Fitted r_B, delta_B and gamma parameters:\n";
  std::cout << "r_B = " << rB << " \u00B1 " << TMath::Sqrt(gammacov(0, 0)) << std::endl;
  std::cout << "delta_B = " << deltaB << " \u00B1 " << TMath::Sqrt(gammacov(1, 1)) << std::endl;
  std::cout << "gamma = " << gamma << " \u00B1 " << TMath::Sqrt(gammacov(2, 2)) << std::endl;
  delete phasespace;
  std::string drawcontours;
  std::cout << "Draw contours?\n";
  std::cin >> drawcontours;
  if(drawcontours == "yes") {
    std::cout << "Drawing contours...\n";
    fitgamma.PlotContours("Contour_rB_vs_dB.png", "Contour_dB_vs_gamma.png", "Contour_gamma_vs_rB.png", 20);
    std::cout << "Finished drawing contours\n";
  }
  std::cout << "Congratulations, gamma has been measured!\n";
  return 0;
}
