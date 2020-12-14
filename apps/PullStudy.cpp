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
#include"TH1D.h"
#include"TMatrixD.h"
#include"KKpipiFit.h"
#include"Constants.h"
#include"BinList.h"
#include"DDecayParameters.h"
#include"CPParameters.h"
#include"Fitter.h"
#include"Gamma.h"
#include"FitGamma.h"
#include"PhaseSpaceParameterisation.h"

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
  TFile *f = new TFile("PullDistributions.root", "RECREATE");
  TTree *PullTree = new TTree("pull", "Pull distributions of #x_#pm and #y_#pm");
  Double_t xplus_pull, xminus_pull, yplus_pull, yminus_pull, rB_pull, dB_pull, gamma_pull, gamma_fitted, gamma_error;
  PullTree->Branch("xplus", &xplus_pull, "xplus/D");
  PullTree->Branch("xminus", &xminus_pull, "xminus/D");
  PullTree->Branch("yplus", &yplus_pull, "yplus/D");
  PullTree->Branch("yminus", &yminus_pull, "yminus/D");
  PullTree->Branch("rB", &rB_pull, "rB/D");
  PullTree->Branch("dB", &dB_pull, "dB/D");
  PullTree->Branch("gamma", &gamma_pull, "gamma/D");
  PullTree->Branch("gamma_fitted", &gamma_fitted, "gamma_fitted/D");
  PullTree->Branch("gamma_error", &gamma_error, "gamma_error/D");
  std::cout << "Tree with pulls is ready\n";
  PhaseSpaceParameterisation *phasespace = KKpipiFit::PickBinningScheme(std::string(argv[1]));
  std::string DDecayFilename = argv[4];
  int SampleSize = atoi(argv[5]);
  int Samples = atoi(argv[6]);
  double gamma_sum = 0, gamma2_sum = 0;
  for(int i = 0; i < Samples; i++) {
    std::cout << "Starting fitting of sample " << i << std::endl;
    BinList binlist(phasespace);
    KKpipiFit::LoadTreesIntoBins(treeBplus, treeBminus, binlist, Samples*i, SampleSize);
    DDecayParameters ddparameters(DDecayFilename);
    Fitter fit(binlist, ddparameters);
    CPParameters cpparameters(0.0, 0.0, 0.0, 0.0);
    double xplus, xminus, yplus, yminus;
    fit.DoFit(cpparameters);
    cpparameters.GetCPParameters(xplus, xminus, yplus, yminus);
    TMatrixD cov = cpparameters.GetCov();
    xplus_pull = (xplus - KKpipi_Constants::xplus)/TMath::Sqrt(cov(0, 0));
    xminus_pull = (xminus - KKpipi_Constants::xminus)/TMath::Sqrt(cov(1, 1));
    yplus_pull = (yplus - KKpipi_Constants::yplus)/TMath::Sqrt(cov(2, 2));
    yminus_pull = (yminus - KKpipi_Constants::yminus)/TMath::Sqrt(cov(3, 3));
    FitGamma gammafitter(cpparameters);
    Gamma gammaparams(0.1, 140, 70);
    gammafitter.DoFit(gammaparams);
    double rB, dB, gamma;
    gammaparams.GetGammaParameters(rB, dB, gamma);
    TMatrixD gammacov = gammaparams.GetCov();
    rB_pull = (rB - KKpipi_Constants::rB)/TMath::Sqrt(gammacov(0, 0));
    dB_pull = (dB - KKpipi_Constants::dB_d)/TMath::Sqrt(gammacov(1, 1));
    gamma_pull = (gamma - KKpipi_Constants::gamma_d)/TMath::Sqrt(gammacov(2, 2));
    gamma_fitted = gamma;
    gamma_error = TMath::Sqrt(gammacov(2, 2));
    PullTree->Fill();
    gamma_sum += gamma;
    gamma2_sum += gamma*gamma;
  }
  std::cout << "Pull study finished\n";
  std::cout << "Analysed " << Samples << " samples, each of size " << SampleSize << std::endl;
  PullTree->Write();
  f->Close();
  fBplus.Close();
  fBminus.Close();
  delete phasespace;
  std::cout << "Standard deviation of gamma: " << TMath::Sqrt((gamma2_sum - gamma_sum*gamma_sum/Samples)/(Samples - 1)) << std::endl;
  std::cout << "Congratulations, gamma has been measured!\n";
  return 0;
}


