// Martin Duy Tat 30th October 2020
/**
 * Fitting is the program for doing the binned fitting
 * D meson decay parameters are calculated using Monte Carlo integration using the ampliude model from AmpGen
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

int main() {
  std::cout << "Starting B->DK, D->KKpipi binned fit\n";
  std::string Bplusfile("/data/lhcb/users/tat/D02KKpipi/BplusEvents/BplustoDK_10K.root");
  std::string Bminusfile("/data/lhcb/users/tat/D02KKpipi/BminusEvents/BminustoDK_10K.root");
  PhaseSpaceParameterisation psp;
  std::cout << "Loaded phase space\n";
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
  double masses[4] = {0.493677, 0.493677, 0.13957039, 0.13957039};
  std::cout << "Generating events for D meson hadronic parameter calculation\n";
  DDecayParameters ddparameters(psp, 1.86483, masses, 1000000);
  std::cout << "Calculated D meson parameters\n";
  Fitter fit(binlist, ddparameters);
  CPParameters cpparameters(-0.09, 0.06, -0.04, 0.08);
  double xplus, xminus, yplus, yminus, xpluse, xminuse, ypluse, yminuse;
  std::cout << "Start fitting\n";
  fit.DoFit(cpparameters);
  std::cout << "Done fitting, getting fitted parameters\n";
  cpparameters.GetCPParameters(xplus, xminus, yplus, yminus);
  cpparameters.GetError(xpluse, xminuse, ypluse, yminuse);
  std::cout << "Fitted parameters:\n";
  std::cout << "xplus = " << xplus << " +- " << xpluse << std::endl;
  std::cout << "xminus = " << xminus << " +- " << xminuse << std::endl;
  std::cout << "yplus = " << yplus << " +- " << ypluse << std::endl;
  std::cout << "yminus = " << yminus << " +- " << yminuse << std::endl;
  return 0;
}
