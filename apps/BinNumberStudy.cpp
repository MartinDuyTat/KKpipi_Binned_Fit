// Martin Duy Tat 30th October 2020
/**
 * Fitting is the program for doing the binned fitting and study how the error depends on the number of bins
 * D meson decay parameters are obtained from an input file
 * @param 1 Filename of B+ event file
 * @param 2 Filename of B- event file
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
#include"RectangularPhaseSpace.h"
#include<stdlib.h>
#include"TCanvas.h"
#include"TGraph.h"
#include"TAxis.h"
#include"TGraphErrors.h"
#include"TLine.h"

int main(int argc, char *argv[]) {
  std::cout << "Starting B->DK, D->KKpipi binned fit\n";
  std::cout << "Using Rectangular binning scheme\n";
  if(argc != 3) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  int N = 9;
  double gamma_unbinned_fit = 74.94, gamma_unbinned_error = 0.32, dB_unbinned_fit = 129.65, dB_unbinned_error = 0.33, rB_unbinned_fit = 0.1011, rB_unbinned_error = 0.0006;
  double xplus_unbinned_fit = rB_unbinned_fit*TMath::Cos((dB_unbinned_fit + gamma_unbinned_fit)*TMath::Pi()/180);
  double xminus_unbinned_fit = rB_unbinned_fit*TMath::Cos((dB_unbinned_fit - gamma_unbinned_fit)*TMath::Pi()/180);
  double yplus_unbinned_fit = rB_unbinned_fit*TMath::Sin((dB_unbinned_fit + gamma_unbinned_fit)*TMath::Pi()/180);
  double yminus_unbinned_fit = rB_unbinned_fit*TMath::Sin((dB_unbinned_fit - gamma_unbinned_fit)*TMath::Pi()/180);
  double xplus_unbinned_error = TMath::Sqrt(TMath::Power(rB_unbinned_error*xplus_unbinned_fit/rB_unbinned_fit, 2)
					    + TMath::Power(dB_unbinned_error*yplus_unbinned_fit*TMath::Pi()/180, 2)
				          + TMath::Power(gamma_unbinned_error*yplus_unbinned_fit*TMath::Pi()/180, 2));
  double xminus_unbinned_error = TMath::Sqrt(TMath::Power(rB_unbinned_error*xminus_unbinned_fit/rB_unbinned_fit, 2)
					   + TMath::Power(dB_unbinned_error*yminus_unbinned_fit*TMath::Pi()/180, 2)
					   + TMath::Power(gamma_unbinned_error*yminus_unbinned_fit*TMath::Pi()/180, 2));
  double yplus_unbinned_error = TMath::Sqrt(TMath::Power(rB_unbinned_error*yplus_unbinned_fit/rB_unbinned_fit, 2)
					  + TMath::Power(dB_unbinned_error*xplus_unbinned_fit*TMath::Pi()/180, 2)
					  + TMath::Power(gamma_unbinned_error*xplus_unbinned_fit*TMath::Pi()/180, 2));
  double yminus_unbinned_error = TMath::Sqrt(TMath::Power(rB_unbinned_error*yminus_unbinned_fit/rB_unbinned_fit, 2)
					   + TMath::Power(dB_unbinned_error*xminus_unbinned_fit*TMath::Pi()/180, 2)
					   + TMath::Power(gamma_unbinned_error*xminus_unbinned_fit*TMath::Pi()/180, 2));
  std::vector<double> NumberBins(N);
  std::vector<double> gamma_fitted(N), gamma_error(N), dB_fitted(N), dB_error(N), rB_fitted(N), rB_error(N);
  std::vector<double> xplus_fitted(N), xplus_error(N), xminus_fitted(N), xminus_error(N), yplus_fitted(N), yplus_error(N), yminus_fitted(N), yminus_error(N);
  std::string Bplusfile = argv[1];
  std::string Bminusfile = argv[2];
  TFile fBplus(Bplusfile.c_str(), "READ");
  TFile fBminus(Bminusfile.c_str(), "READ");
  std::cout << "Opened data files\n";
  TTree *treeBplus, *treeBminus;
  fBplus.GetObject("DalitzEventList", treeBplus);
  fBminus.GetObject("DalitzEventList", treeBminus);
  std::cout << "Loaded trees\n";
  std::cout << "Starting bin number study\n";
  for(int i = 2; i < 11; i++) {
    std::vector<int> bins = {1, 1, 1, i, 2};
    RectangularPhaseSpace phasespace(bins);
    PhaseSpaceParameterisation *psp = &phasespace;
    BinList binlist(psp);
    binlist.LoadTTree(treeBplus, +1);
    binlist.LoadTTree(treeBminus, -1);
    std::string filename = "/data/lhcb/users/tat/D02KKpipi/DHadronicParameters/RectangularPhaseSpace_1_1_1_" + std::to_string(i) + "_2_10M.cs\
v";
    DDecayParameters ddparameters(filename);
    Fitter fit(binlist, ddparameters);
    CPParameters cpparameters(-0.09, 0.06, -0.04, 0.08);
    double xplus, xminus, yplus, yminus;
    fit.DoFit(cpparameters);
    cpparameters.GetCPParameters(xplus, xminus, yplus, yminus);
    TMatrixD cov = cpparameters.GetCov();
    FitGamma fitgamma(cpparameters);
    Gamma gammaparams(0.1, 130.0, 75.0);
    fitgamma.DoFit(gammaparams);
    double rB, deltaB, gamma;
    gammaparams.GetGammaParameters(rB, deltaB, gamma);
    TMatrixD gammacov = gammaparams.GetCov();
    std::cout << "Fitted parameters for " << NumberBins[i - 2] << " bins:\n";
    std::cout << "xplus = " << xplus << " +- " << TMath::Sqrt(cov(0, 0)) << std::endl;
    std::cout << "xminus = " << xminus << " +- " << TMath::Sqrt(cov(1, 1)) << std::endl;
    std::cout << "yplus = " << yplus << " +- " << TMath::Sqrt(cov(2, 2)) << std::endl;
    std::cout << "yminus = " << yminus << " +- " << TMath::Sqrt(cov(3, 3)) << std::endl;
    std::cout << "r_B = " << rB << " +- " << TMath::Sqrt(gammacov(0, 0)) << std::endl;
    std::cout << "delta_B = " << deltaB << " +- " << TMath::Sqrt(gammacov(1, 1)) << std::endl;
    std::cout << "gamma = " << gamma << " +- " << TMath::Sqrt(gammacov(2, 2)) << std::endl;
    NumberBins[i - 2] = i;
    xplus_fitted[i - 2] = xplus;
    xplus_error [i - 2]= TMath::Sqrt(cov(0, 0));
    xminus_fitted[i - 2] = xminus;
    xminus_error[i - 2] = TMath::Sqrt(cov(1, 1));
    yplus_fitted[i - 2] = yplus;
    yplus_error[i - 2] = TMath::Sqrt(cov(2, 2));
    yminus_fitted[i - 2] = yminus;
    yminus_error[i - 2] = TMath::Sqrt(cov(3, 3));
    gamma_fitted[i - 2] = gamma;
    gamma_error[i - 2] = TMath::Sqrt(gammacov(2, 2));
    dB_fitted[i - 2] = deltaB;
    dB_error[i - 2] = TMath::Sqrt(gammacov(1, 1));
    rB_fitted[i - 2] = rB;
    rB_error[i - 2] = TMath::Sqrt(gammacov(0, 0));
  }

  TCanvas *c1 = new TCanvas("gamma", "Number of bins vs #gamma", 200, 10, 600, 400);
  TGraphErrors gr1(N, NumberBins.data(), gamma_fitted.data(), nullptr, gamma_error.data());
  gr1.SetTitle("Number of bins vs #gamma");
  gr1.GetXaxis()->SetTitle("Bins");
  gr1.GetYaxis()->SetTitle("#gamma");
  gr1.Draw("A*");
  TLine gamma_true(NumberBins[0], gamma_unbinned_fit, NumberBins[N - 1], gamma_unbinned_fit);
  gamma_true.SetLineColor(kRed);
  gamma_true.Draw("SAME");
  c1->SaveAs("gamma.png");
  TCanvas *c2 = new TCanvas("gammaerror", "Number of bins vs #gamma error", 200, 10, 600, 400);
  TGraph gr2(N, NumberBins.data(), gamma_error.data());
  gr2.SetTitle("Number of bins vs #gamma error");
  gr2.GetXaxis()->SetTitle("Bins");
  gr2.GetYaxis()->SetTitle("#sigma(#gamma)");
  gr2.Draw("A*");
  TLine gamma_error_true(NumberBins[0], gamma_unbinned_error, NumberBins[N - 1], gamma_unbinned_error);
  gamma_error_true.SetLineColor(kRed);
  gamma_error_true.Draw("SAME");
  gr2.SetMinimum(0.0);
  c2->SaveAs("gammaerror.png");

  TCanvas *c3 = new TCanvas("dB", "Number of bins vs #delta_{B}", 200, 10, 600, 400);
  TGraphErrors gr3(N, NumberBins.data(), dB_fitted.data(), nullptr, dB_error.data());
  gr3.SetTitle("Number of bins vs #delta_{B}");
  gr3.GetXaxis()->SetTitle("Bins");
  gr3.GetYaxis()->SetTitle("#delta_{B}");
  gr3.Draw("A*");
  TLine dB_true(NumberBins[0], dB_unbinned_fit, NumberBins[N - 1], dB_unbinned_fit);
  dB_true.SetLineColor(kRed);
  dB_true.Draw("SAME");
  c3->SaveAs("dB.png");
  TCanvas *c4 = new TCanvas("dBerror", "Number of bins vs #delta_{B} error", 200, 10, 600, 400);
  TGraph gr4(N, NumberBins.data(), dB_error.data());
  gr4.SetTitle("Number of bins vs #delta_{B} error");
  gr4.GetXaxis()->SetTitle("Bins");
  gr4.GetYaxis()->SetTitle("#sigma(#delta_{B})");
  gr4.Draw("A*");
  TLine dB_error_true(NumberBins[0], dB_unbinned_error, NumberBins[N - 1], dB_unbinned_error);
  dB_error_true.SetLineColor(kRed);
  dB_error_true.Draw("SAME");
  gr4.SetMinimum(0.0);
  c4->SaveAs("dBerror.png");

  TCanvas *c5 = new TCanvas("rB", "Number of bins vs r_{B}", 200, 10, 600, 400);
  TGraphErrors gr5(N, NumberBins.data(), rB_fitted.data(), nullptr, rB_error.data());
  gr5.SetTitle("Number of bins vs r_{B}");
  gr5.GetXaxis()->SetTitle("Bins");
  gr5.GetYaxis()->SetTitle("r_{B}");
  gr5.SetMaximum(rB_unbinned_fit*1.05);
  gr5.Draw("A*");
  TLine rB_true(NumberBins[0], rB_unbinned_fit, NumberBins[N - 1], rB_unbinned_fit);
  rB_true.SetLineColor(kRed);
  rB_true.Draw("SAME");
  c5->SaveAs("rB.png");
  TCanvas *c6 = new TCanvas("rBerror", "Number of bins vs r_{B} error", 200, 10, 600, 400);
  TGraph gr6(N, NumberBins.data(), rB_error.data());
  gr6.SetTitle("Number of bins vs r_{B} error");
  gr6.GetXaxis()->SetTitle("Bins");
  gr6.GetYaxis()->SetTitle("#sigma(#r_{B})");
  gr6.Draw("A*");
  TLine rB_error_true(NumberBins[0], rB_unbinned_error, NumberBins[N - 1], rB_unbinned_error);
  rB_error_true.SetLineColor(kRed);
  rB_error_true.Draw("SAME");
  gr6.SetMinimum(0.0);
  c6->SaveAs("rBerror.png");
  
  TCanvas *c7 = new TCanvas("xplus", "Number of bins vs x_{+}", 200, 10, 600, 400);
  TGraphErrors gr7(N, NumberBins.data(), xplus_fitted.data(), nullptr, xplus_error.data());
  gr7.SetTitle("Number of bins vs x_{+}");
  gr7.GetXaxis()->SetTitle("Bins");
  gr7.GetYaxis()->SetTitle("x_{+}");
  gr7.Draw("A*");
  TLine xplus_true(NumberBins[0], xplus_unbinned_fit, NumberBins[N - 1], xplus_unbinned_fit);
  xplus_true.SetLineColor(kRed);
  xplus_true.Draw("SAME");
  c7->SaveAs("xplus.png");
  TCanvas *c8 = new TCanvas("xpluserror", "Number of bins vs x_{+} error", 200, 10, 600, 400);
  TGraph gr8(N, NumberBins.data(), xplus_error.data());
  gr8.SetTitle("Number of bins vs x_{+} error");
  gr8.GetXaxis()->SetTitle("Bins");
  gr8.GetYaxis()->SetTitle("#sigma(x_{+})");
  gr8.Draw("A*");
  TLine xplus_error_true(NumberBins[0], xplus_unbinned_error, NumberBins[N - 1], xplus_unbinned_error);
  xplus_error_true.SetLineColor(kRed);
  xplus_error_true.Draw("SAME");
  gr8.SetMinimum(0.0);
  c8->SaveAs("xpluserror.png");

  TCanvas *c9 = new TCanvas("xminus", "Number of bins vs x_{-}", 200, 10, 600, 400);
  TGraphErrors gr9(N, NumberBins.data(), xminus_fitted.data(), nullptr, xminus_error.data());
  gr9.SetTitle("Number of bins vs x_{-}");
  gr9.GetXaxis()->SetTitle("Bins");
  gr9.GetYaxis()->SetTitle("x_{-}");
  gr9.Draw("A*");
  TLine xminus_true(NumberBins[0], xminus_unbinned_fit, NumberBins[N - 1], xminus_unbinned_fit);
  xminus_true.SetLineColor(kRed);
  xminus_true.Draw("SAME");
  c9->SaveAs("xminus.png");
  TCanvas *c10 = new TCanvas("xminuserror", "Number of bins vs x_{-} error", 200, 10, 600, 400);
  TGraph gr10(N, NumberBins.data(), xminus_error.data());
  gr10.SetTitle("Number of bins vs x_{-} error");
  gr10.GetXaxis()->SetTitle("Bins");
  gr10.GetYaxis()->SetTitle("#sigma(x_{-})");
  gr10.Draw("A*");
  TLine xminus_error_true(NumberBins[0], xminus_unbinned_error, NumberBins[N - 1], xminus_unbinned_error);
  xminus_error_true.SetLineColor(kRed);
  xminus_error_true.Draw("SAME");
  gr10.SetMinimum(0.0);
  c10->SaveAs("xminuserror.png");

  TCanvas *c11 = new TCanvas("yplus", "Number of bins vs y_{+}", 200, 10, 600, 400);
  TGraphErrors gr11(N, NumberBins.data(), yplus_fitted.data(), nullptr, yplus_error.data());
  gr11.SetTitle("Number of bins vs y_{+}");
  gr11.GetXaxis()->SetTitle("Bins");
  gr11.GetYaxis()->SetTitle("y_{+}");
  gr11.Draw("A*");
  TLine yplus_true(NumberBins[0], yplus_unbinned_fit, NumberBins[N - 1], yplus_unbinned_fit);
  yplus_true.SetLineColor(kRed);
  yplus_true.Draw("SAME");
  c11->SaveAs("yplus.png");
  TCanvas *c12 = new TCanvas("ypluserror", "Number of bins vs y_{+} error", 200, 10, 600, 400);
  TGraph gr12(N, NumberBins.data(), yplus_error.data());
  gr12.SetTitle("Number of bins vs y_{+} error");
  gr12.GetXaxis()->SetTitle("Bins");
  gr12.GetYaxis()->SetTitle("#sigma(y_{+})");
  gr12.Draw("A*");
  TLine yplus_error_true(NumberBins[0], yplus_unbinned_error, NumberBins[N - 1], yplus_unbinned_error);
  yplus_error_true.SetLineColor(kRed);
  yplus_error_true.Draw("SAME");
  gr12.SetMinimum(0.0);
  c12->SaveAs("ypluserror.png");

  TCanvas *c13 = new TCanvas("yminus", "Number of bins vs y_{-}", 200, 10, 600, 400);
  TGraphErrors gr13(N, NumberBins.data(), yminus_fitted.data(), nullptr, yminus_error.data());
  gr13.SetTitle("Number of bins vs y_{-}");
  gr13.GetXaxis()->SetTitle("Bins");
  gr13.GetYaxis()->SetTitle("y_{-}");
  gr13.Draw("A*");
  TLine yminus_true(NumberBins[0], yminus_unbinned_fit, NumberBins[N - 1], yminus_unbinned_fit);
  yminus_true.SetLineColor(kRed);
  yminus_true.Draw("SAME");
  c13->SaveAs("yminus.png");
  TCanvas *c14 = new TCanvas("yminuserror", "Number of bins vs y_{-} error", 200, 10, 600, 400);
  TGraph gr14(N, NumberBins.data(), yminus_error.data());
  gr14.SetTitle("Number of bins vs y_{-} error");
  gr14.GetXaxis()->SetTitle("Bins");
  gr14.GetYaxis()->SetTitle("#sigma(y_{-})");
  gr14.Draw("A*");
  TLine yminus_error_true(NumberBins[0], yminus_unbinned_error, NumberBins[N - 1], yminus_unbinned_error);
  yminus_error_true.SetLineColor(kRed);
  yminus_error_true.Draw("SAME");
  gr14.SetMinimum(0.0);
  c14->SaveAs("yminuserror.png");
  return 0;
}
