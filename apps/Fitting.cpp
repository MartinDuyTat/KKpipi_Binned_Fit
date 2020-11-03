// Martin Duy Tat 30th October 2020
/**
 * Fitting is the program for doing the binned fitting
 */

#include<string>
#include<iostream>
#include"PhaseSpaceParameterisation.h"
#include"TTree.h"
#include"BinList.h"
#include"DDecayParameters.h"
#include"CPParameters.h"

int main() {
  std::string Bplusfile("/data/lhcb/users/tat/D02KKpipi/BplusEvents/BplustoDK_1K.root");
  std::string Bminusfile("/data/lhcb/users/tat/D02KKpipi/BminusEvents/BminustoDK_1K.root");
  PhaseSpaceParameterisation psp;
  TFile fBplus(Bplusfile, "READ");
  TFile fBminus(Bminusfile, "READ");
  TTree *treeBplus, *treeBminus;
  fBplus.GetObject("DalitzEventList", treeBplus);
  fBminus.GetObject("DalitzEventList", treeBminus);
  BinList binlist(psp);
  binlist.LoadTree(treeBplus, +1);
  binlist.LoadTree(treeBminus, -1);
  double masses[4] = {0.493677, 0.493677, 0.13957039, 0.13957039};
  DDecayParameters ddparameters(psp, 1.86483, masses, 1000);
  Fitter fit(binlist, ddparameters);
  CPParameters cpparameters(0.1, 0.1, 0.1, 0.1);
  double xplus, xminus, yplus, yminus;
  fit.DoFit(cpparameters);
  cpparameters.GetCPParameters(xplus, xminus, yplus, yminus);
  std::cout << "xplus = " << xplus << std::endl;
  std::cout << "xminus = " << xminus << std::endl;
  std::cout << "yplus = " << yplus << std::endl;
  std::cout << "yminus = " << yminus << std::endl;
  return 0;
}
