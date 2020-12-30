// Martin Duy Tat 19th December 2020
/**
 * BinningQ is the program for calculating the Q-value of a binning scheme
 * D meson decay parameters are obtained from an input file
 * @param 1 Filename of D meson hadronic decay parameters
 */

#include<string>
#include<iostream>
#include"KKpipiMath.h"
#include"TMath.h"

int main(int argc, char *argv[]) {
  if(argc != 2) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  double Q = KKpipiMath::CalculateBinningQValue(std::string(argv[1]));
  double Qexact = KKpipiMath::CalculateExactBinningQValue(std::string(argv[1]));
  std::cout << "Q-value of binning scheme: " << Q << std::endl;
  std::cout << "Q-value using exact formula: " << Qexact << std::endl;
  if(TMath::IsNaN(Q) || Q > 1.0) {
    std::cout << "Binning scheme is invalid!\n";
  } else if(Q < 0.9) {
    std::cout << "Keep working, need to get that Q-value up!\n";
  } else {
    std::cout << "Congratulations, the binning scheme is good!!!\n";
  }
  return 0;
}
