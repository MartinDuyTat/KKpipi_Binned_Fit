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
  BinList binlist(phasespace);
  KKpipiFit::LoadInputDataIntoBins(std::string(argv[2]), std::string(argv[3]), binlist);
  std::cout << "Input events sorted into bins\n";
  std::cout << "Loading hadronic D decay parameters...\n";
  DDecayParameters ddparameters(argv[4]);
  std::cout << "Hadronic D parameters loaded\n";
  std::cout << "Loading fitter...\n";
  Fitter fit(binlist, ddparameters);
  CPParameters cpparameters(0.0, 0.0, 0.0, 0.0);
  std::cout << "Ready for fitting\n";
  std::cout << "Fitting x and y...\n";
  fit.DoFit(cpparameters);
  std::cout << "x and y fitted\n";
  KKpipiFit::PrintXY(cpparameters);
  std::cout << "Fitting r_B, delta_B and gamma...\n";
  FitGamma fitgamma(cpparameters);
  Gamma gammaparams(0.1, 140.0, 70.0);
  fitgamma.DoFit(gammaparams);
  std::cout << "r_B, delta_B and gamma fitted\n";
  KKpipiFit::PrintGamma(gammaparams);
  KKpipiFit::DrawContours(fitgamma);
  std::cout << "Congratulations, gamma has been measured!\n";
  delete phasespace;
  return 0;
}
