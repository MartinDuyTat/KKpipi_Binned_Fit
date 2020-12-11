// Martin Duy Tat 1st December 2020
/**
 * MapStrongPhases is a program that will calculate strong phases on a grid in phase space and save it to a file
 * The loop over i, or x1 bins, has been replaced by a user input to run these jobs in parallel
 * @param 1 Filename for output file with average strong phases
 * @param 2 Index i in the outermost loop over x1 regions
 */

#include<iostream>
#include<fstream>
#include<complex>
#include"KKpipiMath.h"
#include"TRandom3.h"
#include"Constants.h"
#include"Amplitude.h"
#include"SophisticatedPhaseSpace.h"
#include"Event.h"

int main(int argc, char *argv[]) {
  if(argc != 3) {
    std::cout << "Need 2 inputs\n";
    return 0;
  }
  int N = 100;
  std::vector<double> x_low = {KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI, KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI, -1.0, -1.0, -TMath::Pi()};
  std::vector<double> x_high = {KKpipi_Constants::MASS_D - KKpipi_Constants::MASS_K - KKpipi_Constants::MASS_PI, KKpipi_Constants::MASS_D - KKpipi_Constants::MASS_K - KKpipi_Constants::MASS_PI, 1.0, 1.0, TMath::Pi()};
  std::vector<double> dx(5);
  for(int i = 0; i < 5; i++) {
    dx[i] = (x_high[i] - x_low[i])/N;
  }
  Amplitude amplitude;
  SophisticatedPhaseSpace rph(8);
  std::ofstream f(argv[1]);
  int i = atoi(argv[2]);
  if(i >= N) {
    std::cout << "i is outside of range\n";
    return 0;
  }
  if(i == 0) {
    f << N << std::endl;
  }
  for(int j = 0; j < N; j++) {
    for(int k = 0; k < N; k++) {
      std::vector<std::complex<double>> sum(rph.NumberOfRegions(), std::complex<double>(0.0, 0.0));
      for(int n = 0; n < N; n++) {
	for(int m = 0; m < N; m++) {
	  std::vector<double> X = {x_low[0] + dx[0]*(i + 0.5), x_low[1] + dx[1]*(j + 0.5), x_low[2] + dx[2]*(n + 0.5), x_low[3] + dx[3]*(m + 0.5), x_low[4] + dx[4]*(k + 0.5)};
	  std::vector<double> event = KKpipiMath::ConvertXToMomenta(X);
	  int whichregion = rph.WhichRegion(X);
	  //sum[whichregion] += std::arg(amplitude(event, +1))*std::conj(amplitude(event, -1)); //remove maybe
	  sum[whichregion] += std::polar(1.0, std::arg(amplitude(event, +1)/amplitude(event, -1)));
	}
      }
      f << i << "," << j << "," << k << ",";
      for(int bin = 0; bin < rph.NumberOfRegions(); bin++) {
	f << std::arg(sum[bin]) << ",";
      }
      f << std::endl;
    }
  }
  f.close();
  return 0;
}
