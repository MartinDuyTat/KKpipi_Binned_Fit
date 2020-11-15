// Martin Duy Tat 10th November 2020

#include<string>
#include<vector>
#include"FitGamma.h"
#include"Gamma.h"
#include"XYLikelihood.h"
#include"CPParameters.h"
#include"TCanvas.h"
#include"TGraph.h"
#include"TAxis.h"

FitGamma::FitGamma(CPParameters cpparameters): m_cpparameters(cpparameters) {
}

void FitGamma::DoFit(Gamma &GammaParams) {
  if(m_minimiser == nullptr) {
    m_minimiser = new ROOT::Minuit2::Minuit2Minimizer();
  }
  if(m_likelihood == nullptr) {
    m_likelihood = new XYLikelihood(m_cpparameters);
  }
  m_fcn = ROOT::Math::Functor(*m_likelihood, 3);
  m_minimiser->SetFunction(m_fcn);
  double rB, deltaB, gamma;
  GammaParams.GetGammaParameters(rB, deltaB, gamma);
  m_minimiser->SetVariable(0, "r_B", rB, 0.01);
  m_minimiser->SetVariable(1, "delta_B", deltaB, 1);
  m_minimiser->SetVariable(2, "gamma", gamma, 1);
  m_minimiser->Minimize();
  const double *result = m_minimiser->X();
  GammaParams = Gamma(result[0], result[1], result[2]);
  double CovMatrix[16];
  m_minimiser->GetCovMatrix(CovMatrix);
  GammaParams.SetCov(CovMatrix);
}

void FitGamma::PlotContours(std::string Filename_rB_deltaB, std::string Filename_deltaB_gamma, std::string Filename_gamma_rB, unsigned int npoints) const {
  std::vector<double> x(npoints + 1), y(npoints + 1);
  TCanvas *c1 = new TCanvas("rB_vs_dB", "r_{B} vs #delta_{B} contours", 200, 10, 600, 400);
  m_minimiser->SetErrorDef(1.0);
  m_minimiser->Contour(0, 1, npoints, x.data(), y.data());
  x[npoints] = x[0];
  y[npoints] = y[0];
  TGraph *gr1 = new TGraph(npoints + 1, x.data(), y.data());
  gr1->Draw("AC");
  gr1->SetTitle("r_B vs #delta_B contours");
  gr1->GetXaxis()->SetTitle("r_{B}");
  gr1->GetYaxis()->SetTitle("#delta_B");
  c1->SaveAs(Filename_rB_deltaB.c_str());
  delete c1;
  delete gr1;
  TCanvas *c2 = new TCanvas("gamma_vs_deltaB", "#gamma vs #delta_{B} contours", 200, 10, 600, 400);
  m_minimiser->SetErrorDef(1.0);
  m_minimiser->Contour(2, 1, npoints, x.data(), y.data());
  x[npoints] = x[0];
  y[npoints] = y[0];
  TGraph *gr4 = new TGraph(npoints + 1, x.data(), y.data());
  gr4->Draw("AC");
  gr4->SetTitle("#gamma vs #delta_{B} contours");
  gr4->GetXaxis()->SetTitle("#gamma");
  gr4->GetYaxis()->SetTitle("#delta_{B}");
  c2->SaveAs(Filename_deltaB_gamma.c_str());
  delete c2;
  delete gr4;
  TCanvas *c3 = new TCanvas("gamma_vs_rB", "#gamma vs r_{B} contours", 200, 10, 600, 400);
  m_minimiser->SetErrorDef(1.0);
  m_minimiser->Contour(2, 0, npoints, x.data(), y.data());
  x[npoints] = x[0];
  y[npoints] = y[0];
  TGraph *gr7 = new TGraph(npoints + 1, x.data(), y.data());
  gr7->Draw("AC");
  gr7->SetTitle("#gamma vs r_{B} contours");
  gr7->GetXaxis()->SetTitle("#gamma");
  gr7->GetYaxis()->SetTitle("r_{B}");
  c3->SaveAs(Filename_gamma_rB.c_str());
  delete c3;
  delete gr7;
}

FitGamma::~FitGamma() {
  delete m_likelihood;
  m_likelihood = nullptr;
  delete m_minimiser;
  m_minimiser = nullptr;
}
