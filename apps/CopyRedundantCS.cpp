// Martin Duy Tat 5th January 2021
/** 
 * This program is for copying the \f$c_i\f$ and \f$s_i\f$ from a smaller binning scheme to a binning scheme with twice as many bins
 * The parameters are repeated twice because they are constrained to be the same
 * @param 1 Filename to larger binning scheme
 * @param 2 Filename of smaller binning scheme
 * @param 3 Filename of new binning scheme
 */

#include<string>
#include<stdlib.h>
#include<iostream>
#include"DDecayParameters.h"

int main(int argc, char *argv[]) {
  if(argc != 4) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  std::cout << "c_i, s_i transfer from smaller to larger binning scheme\n";
  DDecayParameters larger(argv[1]);
  larger.CopyRedundantCS(std::string(argv[2]));
  larger.SaveCSV(std::string(argv[3]));
  std::cout << "Transfer complete\n";
  return 0;
}
