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

int main(int argc, char *argv[]) {
  std::cout << "Starting B->DK, D->KKpipi binned fit\n";
  if(argc != 4) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  PhaseSpaceParameterisation psp;
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
  std::cout << "Generating events for D meson hadronic parameter calculation\n";
  DDecayParameters ddparameters(argv[3]);
  std::cout << "Calculated D meson parameters\n";
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
  return 0;
}
