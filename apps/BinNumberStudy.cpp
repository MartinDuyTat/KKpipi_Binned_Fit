// Martin Duy Tat 30th October 2020
/**
 * BinNumberStudy is the program for doing the binned fitting and study how the error depends on the number of bins
 * @param 1 Filename of B+ event file with 1M events
 * @param 2 Filename of B- event file with 1M events
 * @param 3 Filename of ROOT file with flat phase space events
 * @param 4 Smallest number of bins
 * @param 5 Largest number of bins
 */

#include<string>
#include<iostream>
#include<stdlib.h>
#include<numeric>
#include"TFile.h"
#include"TTree.h"
#include"TMath.h"
#include"TMatrixD.h"
#include"TCanvas.h"
#include"TGraph.h"
#include"TAxis.h"
#include"TGraphErrors.h"
#include"TMultiGraph.h"
#include"TLine.h"
#include"KKpipiFit.h"
#include"PhaseSpaceParameterisation.h"
#include"AmplitudePhaseSpace.h"
#include"BinList.h"
#include"DDecayParameters.h"
#include"CPParameters.h"
#include"Fitter.h"
#include"FitGamma.h"
#include"Gamma.h"
#include"EventList.h"

// Results from unbinned fit with 1M events of each sample
const static double gamma_unbinned_fit   = 74.94,
                    gamma_unbinned_error = 0.32,
                    dB_unbinned_fit      = 129.65,
                    dB_unbinned_error    = 0.33,
                    rB_unbinned_fit      = 0.1011,
                    rB_unbinned_error    = 0.0006;

// Use definition of x and y
const static double xplus_unbinned_fit  = rB_unbinned_fit*TMath::Cos((dB_unbinned_fit + gamma_unbinned_fit)*TMath::Pi()/180),
                    xminus_unbinned_fit = rB_unbinned_fit*TMath::Cos((dB_unbinned_fit - gamma_unbinned_fit)*TMath::Pi()/180),
                    yplus_unbinned_fit  = rB_unbinned_fit*TMath::Sin((dB_unbinned_fit + gamma_unbinned_fit)*TMath::Pi()/180),
                    yminus_unbinned_fit = rB_unbinned_fit*TMath::Sin((dB_unbinned_fit - gamma_unbinned_fit)*TMath::Pi()/180);

// Error propagation, assuming errors are uncorrelated (alright approximation)
const static double xplus_unbinned_error  = TMath::Sqrt(TMath::Power(rB_unbinned_error*xplus_unbinned_fit/rB_unbinned_fit, 2)
					              + TMath::Power(dB_unbinned_error*yplus_unbinned_fit*TMath::Pi()/180, 2)
				                      + TMath::Power(gamma_unbinned_error*yplus_unbinned_fit*TMath::Pi()/180, 2));

const static double xminus_unbinned_error = TMath::Sqrt(TMath::Power(rB_unbinned_error*xminus_unbinned_fit/rB_unbinned_fit, 2)
					              + TMath::Power(dB_unbinned_error*yminus_unbinned_fit*TMath::Pi()/180, 2)
					              + TMath::Power(gamma_unbinned_error*yminus_unbinned_fit*TMath::Pi()/180, 2));

const static double yplus_unbinned_error  = TMath::Sqrt(TMath::Power(rB_unbinned_error*yplus_unbinned_fit/rB_unbinned_fit, 2)
					              + TMath::Power(dB_unbinned_error*xplus_unbinned_fit*TMath::Pi()/180, 2)
              					      + TMath::Power(gamma_unbinned_error*xplus_unbinned_fit*TMath::Pi()/180, 2));

const static double yminus_unbinned_error = TMath::Sqrt(TMath::Power(rB_unbinned_error*yminus_unbinned_fit/rB_unbinned_fit, 2)
					              + TMath::Power(dB_unbinned_error*xminus_unbinned_fit*TMath::Pi()/180, 2)
					              + TMath::Power(gamma_unbinned_error*xminus_unbinned_fit*TMath::Pi()/180, 2));

int main(int argc, char *argv[]) {
  std::cout << "Starting B->DK, D->KKpipi binned fit bin number study\n";
  std::cout << "Using Amplitude binning scheme\n";
  if(argc != 6) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  int StartBin = atoi(argv[4]);
  int EndBin = atoi(argv[5]);
  int N = EndBin - StartBin + 1;
  std::vector<double> NumberBins(N);
  std::vector<double> gamma_fitted, gamma_error, dB_fitted, dB_error, rB_fitted, rB_error;
  std::vector<double> xplus_fitted, xplus_error, xminus_fitted, xminus_error, yplus_fitted, yplus_error, yminus_fitted, yminus_error;
  std::cout << "Loading input events...\n";
  TFile fBplus(argv[1], "READ");
  TFile fBminus(argv[2], "READ");
  std::cout << "Opened event files\n";
  std::cout << "Loading events into trees...\n";
  TTree *treeBplus, *treeBminus;
  fBplus.GetObject("DalitzEventList", treeBplus);
  fBminus.GetObject("DalitzEventList", treeBminus);
  std::cout << "Loaded trees\n";
  std::cout << "Loading flat phase space events...\n";
  EventList eventlist;
  eventlist.LoadTree(std::string(argv[3]));
  std::cout << "Flat phase space events ready\n";
  std::iota(NumberBins.begin(), NumberBins.end(), StartBin);
  std::cout << "Starting bin number study\n";
  for(const auto &i : NumberBins) {
    AmplitudePhaseSpace phasespace(i);
    PhaseSpaceParameterisation *psp = &phasespace;
    BinList binlist(psp);
    KKpipiFit::LoadTreesIntoBins(treeBplus, treeBminus, binlist);
    phasespace.ReadAmplitudeFromEvent(true);
    DDecayParameters ddparameters(psp, eventlist);
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
    std::cout << "Fitted parameters for " << i << " bins:\n";
    std::cout << "xplus = " << xplus << " \u00b1 " << TMath::Sqrt(cov(0, 0)) << std::endl;
    std::cout << "xminus = " << xminus << " \u00b1 " << TMath::Sqrt(cov(1, 1)) << std::endl;
    std::cout << "yplus = " << yplus << " \u00b1 " << TMath::Sqrt(cov(2, 2)) << std::endl;
    std::cout << "yminus = " << yminus << " \u00b1 " << TMath::Sqrt(cov(3, 3)) << std::endl;
    std::cout << "r_B = " << rB << " \u00b1 " << TMath::Sqrt(gammacov(0, 0)) << std::endl;
    std::cout << "delta_B = " << deltaB << " \u00b1 " << TMath::Sqrt(gammacov(1, 1)) << std::endl;
    std::cout << "gamma = " << gamma << " \u00b1 " << TMath::Sqrt(gammacov(2, 2)) << std::endl;
    xplus_fitted.push_back(xplus);
    xplus_error.push_back(TMath::Sqrt(cov(0, 0)));
    xminus_fitted.push_back(xminus);
    xminus_error.push_back(TMath::Sqrt(cov(1, 1)));
    yplus_fitted.push_back(yplus);
    yplus_error.push_back(TMath::Sqrt(cov(2, 2)));
    yminus_fitted.push_back(yminus);
    yminus_error.push_back(TMath::Sqrt(cov(3, 3)));
    gamma_fitted.push_back(gamma);
    gamma_error.push_back(TMath::Sqrt(gammacov(2, 2)));
    dB_fitted.push_back(deltaB);
    dB_error.push_back(TMath::Sqrt(gammacov(1, 1)));
    rB_fitted.push_back(rB);
    rB_error.push_back(TMath::Sqrt(gammacov(0, 0)));
  }

  TCanvas c1("gamma", "Number of bins vs #gamma", 200, 10, 600, 400);
  TMultiGraph mg1("mg1", "mg1");
  mg1.SetTitle("Number of bins vs #gamma;Bins;#gamma");
  TGraphErrors gr1(N, NumberBins.data(), gamma_fitted.data(), nullptr, gamma_error.data());
  gr1.SetTitle("Binned fit of #gamma");
  mg1.Add(&gr1, "*");
  double x1line[2] = {NumberBins[0], NumberBins.back()}, y1line[2] = {gamma_unbinned_fit, gamma_unbinned_fit};
  TGraph unbinned1(2, x1line, y1line);
  unbinned1.SetLineColor(kRed);
  unbinned1.SetTitle("Unbinned fit of #gamma");
  mg1.Add(&unbinned1, "L");
  mg1.Draw("A");
  c1.SaveAs("gamma.png");

  TCanvas c2("gammaerror", "Number of bins vs #gamma error", 200, 10, 600, 400);
  TMultiGraph mg2("mg2", "mg2");
  mg2.SetTitle("Number of bins vs #gamma error;Bins;#sigma(#gamma)");
  TGraph gr2(N, NumberBins.data(), gamma_error.data());
  gr2.SetTitle("Binned fit error of #gamma");
  mg2.Add(&gr2, "*");
  double x2line[2] = {NumberBins[0], NumberBins.back()}, y2line[2] = {gamma_unbinned_error, gamma_unbinned_error};
  TGraph unbinned2(2, x2line, y2line);
  unbinned2.SetLineColor(kRed);
  unbinned2.SetTitle("Unbinned fit error of #gamma");
  mg2.Add(&unbinned2, "L");
  mg2.SetMinimum(0.0);
  mg2.Draw("A");
  c2.SaveAs("gammaerror.png");

  TCanvas c3("dB", "Number of bins vs #delta_{B}", 200, 10, 600, 400);
  TMultiGraph mg3("mg3", "mg3");
  mg3.SetTitle("Number of bins vs #delta_{B};Bins;#delta_{B}");
  TGraphErrors gr3(N, NumberBins.data(), dB_fitted.data(), nullptr, dB_error.data());
  gr3.SetTitle("Binned fit of #delta_{B}");
  mg3.Add(&gr3, "*");
  double x3line[2] = {NumberBins[0], NumberBins.back()}, y3line[2] = {dB_unbinned_fit, dB_unbinned_fit};
  TGraph unbinned3(2, x3line, y3line);
  unbinned3.SetLineColor(kRed);
  unbinned3.SetTitle("Unbinned fit of #delta_{B}");
  mg3.Add(&unbinned3, "L");
  mg3.Draw("A");
  c3.SaveAs("dB.png");

  TCanvas c4("dBerror", "Number of bins vs #delta_{B} error", 200, 10, 600, 400);
  TMultiGraph mg4("mg4", "mg4");
  mg4.SetTitle("Number of bins vs #delta_{B} error;Bins;#sigma(#delta_{B})");
  TGraph gr4(N, NumberBins.data(), dB_error.data());
  gr4.SetTitle("Binned fit error of #delta_{B}");
  mg4.Add(&gr4, "*");
  double x4line[2] = {NumberBins[0], NumberBins.back()}, y4line[2] = {dB_unbinned_error, dB_unbinned_error};
  TGraph unbinned4(2, x4line, y4line);
  unbinned4.SetLineColor(kRed);
  unbinned4.SetTitle("Unbinned fit error of #delta_{B}");
  mg4.Add(&unbinned4, "L");
  mg4.SetMinimum(0.0);
  mg4.Draw("A");
  c4.SaveAs("dBerror.png");

  TCanvas c5("rB", "Number of bins vs r_{B}", 200, 10, 600, 400);
  TMultiGraph mg5("mg5", "mg5");
  mg5.SetTitle("Number of bins vs r_{B};Bins;r_{B}");
  TGraphErrors gr5(N, NumberBins.data(), rB_fitted.data(), nullptr, rB_error.data());
  gr5.SetTitle("Binned fit of r_{B}");
  mg5.Add(&gr5, "*");
  double x5line[2] = {NumberBins[0], NumberBins.back()}, y5line[2] = {rB_unbinned_fit, rB_unbinned_fit};
  TGraph unbinned5(2, x5line, y5line);
  unbinned5.SetLineColor(kRed);
  unbinned5.SetTitle("Unbinned fit of r_{B}");
  mg5.Add(&unbinned5, "L");
  mg5.Draw("A");
  c5.SaveAs("rB.png");

  TCanvas c6("rBerror", "Number of bins vs r_{B} error", 200, 10, 600, 400);
  TMultiGraph mg6("mg6", "mg6");
  mg6.SetTitle("Number of bins vs r_{B} error;Bins;#sigma(r_{B})");
  TGraph gr6(N, NumberBins.data(), rB_error.data());
  gr6.SetTitle("Binned fit error of r_{B}");
  mg6.Add(&gr6, "*");
  double x6line[2] = {NumberBins[0], NumberBins.back()}, y6line[2] = {rB_unbinned_error, rB_unbinned_error};
  TGraph unbinned6(2, x6line, y6line);
  unbinned6.SetLineColor(kRed);
  unbinned6.SetTitle("Unbinned fit error of r_{B}");
  mg6.Add(&unbinned6, "L");
  mg6.SetMinimum(0.0);
  mg6.Draw("A");
  c6.SaveAs("rBerror.png");
  
  TCanvas c7("xplus", "Number of bins vs x_{+}", 200, 10, 600, 400);
  TMultiGraph mg7("mg7", "mg7");
  mg7.SetTitle("Number of bins vs x_{+};Bins;x_{+}");
  TGraphErrors gr7(N, NumberBins.data(), xplus_fitted.data(), nullptr, xplus_error.data());
  gr7.SetTitle("Binned fit of x_{+}");
  mg7.Add(&gr7, "*");
  double x7line[2] = {NumberBins[0], NumberBins.back()}, y7line[2] = {xplus_unbinned_fit, xplus_unbinned_fit};
  TGraph unbinned7(2, x7line, y7line);
  unbinned7.SetLineColor(kRed);
  unbinned7.SetTitle("Unbinned fit of x_{+}");
  mg7.Add(&unbinned7, "L");
  mg7.Draw("A");
  c7.SaveAs("xplus.png");

  TCanvas c8("xpluserror", "Number of bins vs x_{+} error", 200, 10, 600, 400);
  TMultiGraph mg8("mg8", "mg8");
  mg8.SetTitle("Number of bins vs x_{+} error;Bins;#sigma(x_{+})");
  TGraph gr8(N, NumberBins.data(), xplus_error.data());
  gr8.SetTitle("Binned fit error of x_{+}");
  mg8.Add(&gr8, "*");
  double x8line[2] = {NumberBins[0], NumberBins.back()}, y8line[2] = {xplus_unbinned_error, xplus_unbinned_error};
  TGraph unbinned8(2, x8line, y8line);
  unbinned8.SetLineColor(kRed);
  unbinned8.SetTitle("Unbinned fit error of x_{+}");
  mg8.Add(&unbinned8, "L");
  mg8.SetMinimum(0.0);
  mg8.Draw("A");
  c8.SaveAs("xpluserror.png");

  TCanvas c9("xminus", "Number of bins vs x_{-}", 200, 10, 600, 400);
  TMultiGraph mg9("mg9", "mg9");
  mg9.SetTitle("Number of bins vs x_{-};Bins;x_{-}");
  TGraphErrors gr9(N, NumberBins.data(), xminus_fitted.data(), nullptr, xminus_error.data());
  gr9.SetTitle("Binned fit of x_{-}");
  mg9.Add(&gr9, "*");
  double x9line[2] = {NumberBins[0], NumberBins.back()}, y9line[2] = {xminus_unbinned_fit, xminus_unbinned_fit};
  TGraph unbinned9(2, x9line, y9line);
  unbinned9.SetLineColor(kRed);
  unbinned9.SetTitle("Unbinned fit of x_{-}");
  mg9.Add(&unbinned9, "L");
  mg9.Draw("A");
  c9.SaveAs("xminus.png");

  TCanvas c10("xminuserror", "Number of bins vs x_{-} error", 200, 10, 600, 400);
  TMultiGraph mg10("mg10", "mg10");
  mg10.SetTitle("Number of bins vs x_{-} error;Bins;#sigma(x_{-})");
  TGraph gr10(N, NumberBins.data(), xminus_error.data());
  gr10.SetTitle("Binned fit error of x_{-}");
  mg10.Add(&gr10, "*");
  double x10line[2] = {NumberBins[0], NumberBins.back()}, y10line[2] = {xminus_unbinned_error, xminus_unbinned_error};
  TGraph unbinned10(2, x10line, y10line);
  unbinned10.SetLineColor(kRed);
  unbinned10.SetTitle("Unbinned fit error of x_{-}");
  mg10.Add(&unbinned10, "L");
  mg10.SetMinimum(0.0);
  mg10.Draw("A");
  c10.SaveAs("xminuserror.png");

  TCanvas c11("yplus", "Number of bins vs y_{+}", 200, 10, 600, 400);
  TMultiGraph mg11("mg11", "mg11");
  mg11.SetTitle("Number of bins vs y_{+};Bins;y_{+}");
  TGraphErrors gr11(N, NumberBins.data(), yplus_fitted.data(), nullptr, yplus_error.data());
  gr11.SetTitle("Binned Fit of y_{+}");
  mg11.Add(&gr11, "*");
  double x11line[2] = {NumberBins[0], NumberBins.back()}, y11line[2] = {yplus_unbinned_fit, yplus_unbinned_fit};
  TGraph unbinned11(2, x11line, y11line);
  unbinned11.SetLineColor(kRed);
  unbinned11.SetTitle("Unbinned fit of y_{+}");
  mg11.Add(&unbinned11, "L");
  mg11.Draw("A");
  c11.SaveAs("yplus.png");

  TCanvas c12("ypluserror", "Number of bins vs y_{+} error", 200, 10, 600, 400);
  TMultiGraph mg12("mg12", "mg12");
  mg12.SetTitle("Number of bins vs y_{+} error;Bins;#sigma(y_{+})");
  TGraph gr12(N, NumberBins.data(), yplus_error.data());
  gr12.SetTitle("Binned fit error of y_{+}");
  mg12.Add(&gr12, "*");
  double x12line[2] = {NumberBins[0], NumberBins.back()}, y12line[2] = {yplus_unbinned_error, yplus_unbinned_error};
  TGraph unbinned12(2, x12line, y12line);
  unbinned12.SetLineColor(kRed);
  unbinned12.SetTitle("Unbinned fit error of y_{+}");
  mg12.Add(&unbinned12, "L");
  mg12.SetMinimum(0.0);
  mg12.Draw("A");
  c12.SaveAs("ypluserror.png");

  TCanvas c13("yminus", "Number of bins vs y_{-}", 200, 10, 600, 400);
  TMultiGraph mg13("mg13", "mg13");
  mg13.SetTitle("Number of bins vs y_{-};Bins;y_{-}");
  TGraphErrors gr13(N, NumberBins.data(), yminus_fitted.data(), nullptr, yminus_error.data());
  gr13.SetTitle("Binned fit of y_{-}");
  mg13.Add(&gr13, "*");
  double x13line[2] = {NumberBins[0], NumberBins.back()}, y13line[2] = {yminus_unbinned_fit, yminus_unbinned_fit};
  TGraph unbinned13(2, x13line, y13line);
  unbinned13.SetLineColor(kRed);
  unbinned13.SetTitle("Unbinned fit of y_{-}");
  mg13.Add(&unbinned13, "L");
  mg13.Draw("A");
  c13.SaveAs("yminus.png");

  TCanvas c14("yminuserror", "Number of bins vs y_{-} error", 200, 10, 600, 400);
  TMultiGraph mg14("mg14", "mg14");
  mg14.SetTitle("Number of bins vs y_{-} error;Bins;#sigma(y_{-})");
  TGraph gr14(N, NumberBins.data(), yminus_error.data());
  gr14.SetTitle("Binned fit error of y_{-}");
  mg14.Add(&gr14, "*");
  double x14line[2] = {NumberBins[0], NumberBins.back()}, y14line[2] = {yminus_unbinned_error, yminus_unbinned_error};
  TGraph unbinned14(2, x14line, y14line);
  unbinned14.SetLineColor(kRed);
  unbinned14.SetTitle("Unbinned fit error of y_{-}");
  mg14.Add(&unbinned14, "L");
  mg14.SetMinimum(0.0);
  mg14.Draw("A");
  c14.SaveAs("yminuserror.png");

  return 0;
}
