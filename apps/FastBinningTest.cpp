// Martin Duy Tat 30th October 2020
/**
 * Pullstudy is the program for doing a pull study of the binned fitting
 * D meson decay parameters are loaded from file
 * @param 1 Binning choice, in lowercase letters
 * @param 2 Filename of B+ event file
 * @param 3 Filename of B- event file
 * @param 4 Filename of D meson hadronic decay parameters
 * @param 5 Sample size
 * @param 6 Number of samples
 */

#include<string>
#include<iostream>
#include<stdlib.h>
#include"TFile.h"
#include"TTree.h"
#include"TMath.h"
#include"TMatrixD.h"
#include"KKpipiFit.h"
#include"BinList.h"
#include"DDecayParameters.h"
#include"CPParameters.h"
#include"Fitter.h"
#include"Gamma.h"
#include"FitGamma.h"
#include"PhaseSpaceParameterisation.h"
//#include"SophisticatedPhaseSpace.h"
//#include"AmplitudePhaseSpace.h"
//#include"NaivePhaseSpace.h"

int main(int argc, char *argv[]) {
  if(argc != 7) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  std::cout << "Starting B->DK, D->KKpipi binned fit pull study\n";
  std::cout << "Loading input data...\n";
  TFile fBplus(argv[2], "READ");
  TFile fBminus(argv[3], "READ");
  TTree *treeBplus = nullptr, *treeBminus = nullptr;
  fBplus.GetObject("DalitzEventList", treeBplus);
  fBminus.GetObject("DalitzEventList", treeBminus);
  std::cout << "Events loaded into bins\n";
  std::cout << "Loading pull tree...\n";
  std::cout << "Tree with pulls is ready\n";
  PhaseSpaceParameterisation *phasespace = KKpipiFit::PickBinningScheme(std::string(argv[1]));
  std::string DDecayFilename = argv[4];
  int SampleSize = atoi(argv[5]);
  int Samples = atoi(argv[6]);
  double gamma_sum = 0, gamma2_sum = 0;
  for(int i = 0; i < Samples; i++) {
    BinList binlist(phasespace);
    KKpipiFit::LoadTreesIntoBins(treeBplus, treeBminus, binlist, Samples*i, SampleSize);
    DDecayParameters ddparameters(DDecayFilename);
    Fitter fit(binlist, ddparameters);
    CPParameters cpparameters(0.0, 0.0, 0.0, 0.0);
    fit.DoFit(cpparameters);
    FitGamma gammafitter(cpparameters);
    Gamma gammaparams(0.1, 140, 70);
    gammafitter.DoFit(gammaparams);
    double rB, dB, gamma;
    gammaparams.GetGammaParameters(rB, dB, gamma);
    TMatrixD gammacov = gammaparams.GetCov();
    gamma_sum += gamma;
    gamma2_sum += gamma*gamma;
  }
  std::cout << "Pull study finished\n";
  std::cout << "Analysed " << Samples << " samples, each of size " << SampleSize << std::endl;
  delete phasespace;
  std::cout << "Mean of gamma:               " << gamma_sum/Samples << std::endl;
  std::cout << "Standard deviation of gamma: " << TMath::Sqrt((gamma2_sum - gamma_sum*gamma_sum/Samples)/(Samples - 1)) << std::endl;
  std::cout << "Congratulations, gamma has been measured!\n";
  return 0;
}


