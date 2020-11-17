// Martin Duy Tat 30th October 2020
/**
 * Fitting is the program for doing the binned fitting
 * D meson decay parameters are obtained from an input file
 * @param 1 Filename of B+ event file
 * @param 2 Filename of B- event file
 * @param 3 Filename of D meson hadronic decay parameters
 */

#include<string>
#include<iostream>
#include"PhaseSpaceParameterisation.h"
#include"TFile.h"
#include"TTree.h"
#include"BinList.h"
#include"DDecayParameters.h"
#include"CPParameters.h"
#include"Fitter.h"
#include"TMath.h"
#include"FitGamma.h"
#include"Gamma.h"
#include"TMatrixD.h"
#include"NaiivePhaseSpace.h"
#include"RectangularPhaseSpace.h"

int main(int argc, char *argv[]) {
  std::cout << "Starting B->DK, D->KKpipi binned fit\n";
  if(argc != 4) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  std::vector<int> bins = {2, 2, 2, 2, 2};
  RectangularPhaseSpace phasespace(bins);
  PhaseSpaceParameterisation *psp = &phasespace;
  std::cout << "Loaded phase space\n";
  std::string Bplusfile = argv[1];
  std::string Bminusfile = argv[2];
  TFile fBplus(Bplusfile.c_str(), "READ");
  TFile fBminus(Bminusfile.c_str(), "READ");
  std::cout << "Opened data files\n";
  TTree *treeBplus, *treeBminus;
  fBplus.GetObject("DalitzEventList", treeBplus);
  fBminus.GetObject("DalitzEventList", treeBminus);
  std::cout << "Loaded trees\n";
  BinList binlist(psp);
  std::cout << "Loaded bins\n";
  binlist.LoadTTree(treeBplus, +1);
  binlist.LoadTTree(treeBminus, -1);
  std::cout << "Loaded tree events into bins\n";
  DDecayParameters ddparameters(argv[3]);
  std::cout << "Loaded D meson hadronic parameters\n";
  Fitter fit(binlist, ddparameters);
  CPParameters cpparameters(-0.09, 0.06, -0.04, 0.08);
  double xplus, xminus, yplus, yminus;
  std::cout << "Start fitting\n";
  fit.DoFit(cpparameters);
  std::cout << "Done fitting, getting fitted parameters\n";
  cpparameters.GetCPParameters(xplus, xminus, yplus, yminus);
  TMatrixD cov = cpparameters.GetCov();
  std::cout << "Fitted parameters:\n";
  std::cout << "xplus = " << xplus << " +- " << TMath::Sqrt(cov(0, 0)) << std::endl;
  std::cout << "xminus = " << xminus << " +- " << TMath::Sqrt(cov(1, 1)) << std::endl;
  std::cout << "yplus = " << yplus << " +- " << TMath::Sqrt(cov(2, 2)) << std::endl;
  std::cout << "yminus = " << yminus << " +- " << TMath::Sqrt(cov(3, 3)) << std::endl;
  std::cout << "Starting fit to determine r_B, delta_B and gamma\n";
  FitGamma fitgamma(cpparameters);
  Gamma gammaparams(0.05, 140.0, 60.0);
  fitgamma.DoFit(gammaparams);
  std::cout << "Done fitting for r_B, delta_B and gamma\n";
  double rB, deltaB, gamma;
  gammaparams.GetGammaParameters(rB, deltaB, gamma);
  TMatrixD gammacov = gammaparams.GetCov();
  std::cout << "Fitted parameters:\n";
  std::cout << "r_B = " << rB << " +- " << TMath::Sqrt(gammacov(0, 0)) << std::endl;
  std::cout << "delta_B = " << deltaB << " +- " << TMath::Sqrt(gammacov(1, 1)) << std::endl;
  std::cout << "gamma = " << gamma << " +- " << TMath::Sqrt(gammacov(2, 2)) << std::endl;
  std::cout << "Drawing contours\n";
  fitgamma.PlotContours("Contour_rB_vs_dB.png", "Contour_dB_vs_gamma.png", "Contour_gamma_vs_rB.png", 20);
  std::cout << "Finished drawing contours\n";
  return 0;
}
