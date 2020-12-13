// Martin Duy Tat 30th October 2020
/**
 * Pullstudy is the program for doing a pull study of the binned fitting
 * D meson decay parameters are loaded from file
 * @param 1 Filename of B+ event file
 * @param 2 Filename of B- event file
 * @param 3 Filename of D meson hadronic decay parameters
 * @param 4 Sample size
 * @param 5 Number of samples
 * @param 6 Filename of mean phases in the \f$(x_1, x_2, x_5)\f$ volume
 * @param 7 Total number of bins
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
#include"SophisticatedPhaseSpace.h"
#include"AmplitudePhaseSpace.h"
#include"NaivePhaseSpace.h"

void SplitTree(TTree *tree, TTree *treeSmall, const int &StartEvent, const int &SampleSize);

int main(int argc, char *argv[]) {
  if(argc != 7) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  std::string BplusFilename = argv[1];
  std::string BminusFilename = argv[2];
  std::string DDecayFilename = argv[3];
  int SampleSize = atoi(argv[4]);
  int Samples = atoi(argv[5]);
  std::cout << "Starting B->DK, D->KKpipi binned fit pool study\n";
  std::cout << "Loaded phase space\n";
  TFile fBplus(BplusFilename.c_str(), "READ");
  TFile fBminus(BminusFilename.c_str(), "READ");
  std::cout << "Opened data files\n";
  TTree *treeBplus, *treeBminus;
  fBplus.GetObject("DalitzEventList", treeBplus);
  fBminus.GetObject("DalitzEventList", treeBminus);
  std::cout << "Loaded trees\n";
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
  //SophisticatedPhaseSpace phasespace(atoi(argv[7]));
  AmplitudePhaseSpace phasespace(atoi(argv[6]));
  //phasespace.ReadAverageStrongPhases(std::string(argv[6]));
  PhaseSpaceParameterisation *psp = &phasespace;
  for(int i = 0; i < Samples; i++) {
    std::cout << "Starting fitting of sample " << i << std::endl;
    TTree *treeSmallBplus = new TTree("DalitzEventList", "Dbar0 K+ K- pi+ pi-");
    TTree *treeSmallBminus = new TTree("DalitzEventList", "D0 K+ K- pi+ pi-");
    KKpipiFit::SplitTree(treeBplus, treeSmallBplus, Samples*i, SampleSize);
    KKpipiFit::SplitTree(treeBminus, treeSmallBminus, Samples*i, SampleSize);
    BinList binlist(psp);
    binlist.LoadTTree(treeSmallBplus, +1);
    binlist.LoadTTree(treeSmallBminus, -1);    
    DDecayParameters ddparameters(DDecayFilename);
    Fitter fit(binlist, ddparameters);
    //CPParameters cpparameters(-0.09, 0.06, -0.04, 0.08);
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
    //delete treeSmallBplus;
    //delete treeSmallBminus;
  }
  std::cout << "Pull study finished\n";
  std::cout << "Analysed " << Samples << " samples, each of size " << SampleSize << std::endl;
  PullTree->Write();
  f->Close();
  return 0;
}


