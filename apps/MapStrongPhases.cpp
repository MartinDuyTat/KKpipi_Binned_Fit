// Martin Duy Tat 1st December 2020
/**
 * MapStrongPhases is a program that will calculate strong phases on a grid in phase space and save it to a file
 */

#include<iostream>
#include<fstream>
#include<complex>
#include"KKpipiMath.h"
#include"TRandom3.h"
#include"Constants.h"
#include"Amplitude.h"

int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cout << "Need 5 inputs\n";
    return 0;
  }
  int N = 20;
  std::vector<double> x_low = {KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI, KKpipi_Constants::MASS_K + KKpipi_Constants::MASS_PI, -1.0, -1.0, -TMath::Pi()};
  std::vector<double> x_high = {KKpipi_Constants::MASS_D - KKpipi_Constants::MASS_K - KKpipi_Constants::MASS_PI, KKpipi_Constants::MASS_D - KKpipi_Constants::MASS_K - KKpipi_Constants::MASS_PI, 1.0, 1.0, TMath::Pi()};
  std::vector<double> dx(5);
  for(int i = 0; i < 5; i++) {
    dx[i] = (x_high[i] - x_low[i])/N;
  }
  Amplitude amplitude("D0toKKpipi.so", "Dbar0toKKpipi.so");
  std::ofstream f(argv[1]);
  f << N << std::endl;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      for(int k = 0; k < N; k++) {
	std::complex<double> sum(0.0, 0.0);
	int M = 0;
	for(int n = 0; n < N; n++) {
	  for(int m = 0; m < N; m++) {
	    std::vector<double> X = {x_low[0] + dx[0]*(i + 0.5), x_low[1] + dx[1]*(j + 0.5), x_low[2] + dx[2]*(n + 0.5), x_low[3] + dx[3]*(m + 0.5), x_low[4] + dx[4]*(k + 0.5)};
	    if(X[2] + X[3] < 1.4) {
	      continue;
	    }
	    ++M;
	    std::vector<double> event = KKpipiMath::ConvertXToMomenta(X);
	    sum += std::polar(1.0, std::arg(amplitude(event, +1)/amplitude(event, -1)));
	  }
	}
	f << i << "," << j << "," << k << "," << std::arg(sum/(M*1.0)) << std::endl;
      }
    }
  }
  f.close();
  return 0;
}
