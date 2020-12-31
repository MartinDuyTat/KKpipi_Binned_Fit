// Martin Duy Tat 31st December 2020
/**
 * BinningOptimizer optimizes the binning \f$Q\f$-value
 * D hadronic parameters are calculated from flat phase space events
 * @param 1 Number of bins, must be an even number
 * @param 2 Filename of ROOT file with flat phase space events
 */

#include<string>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include"AmplitudeBinningOptimizer.h"

int main(int argc, char *argv[]) {
  if(argc != 3) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  std::cout << "Optimizing bin edges for Amplitude binning scheme\n";
  std::cout << "Loading flat phase space events...\n";
  AmplitudeBinningOptimizer optimizer(atoi(argv[1]), std::string(argv[2]));
  std::cout << "Events loaded into optimizer\n";
  std::cout << "Optimizing...\n";
  double OptimizedQ = optimizer.OptimizeBinEdges();
  std::cout << "Bin edges optimized\n";
  std::cout << "Optimized binning Q-value: " << OptimizedQ << std::endl;
  std::vector<double> BinEdges = optimizer.GetBinEdges();
  std::cout << "Optimized bin edges: ";
  for(const auto &edge : BinEdges) {
    std::cout << edge << " ";
  }
  std::cout << std::endl;
  return 0;
}
