// Martin Duy Tat 9th November 2020
/** 
 * This program is for plotting hadronic parameters of D->KKpipi decays, and the parameters are read from a CSV file
 * @param 1 Filename of D hadronic decay parameters
 * @param 2 Filename of s_i vs c_i plot
 * @param 3 Filename of K_i and Kbar_i plot
 */

#include<string>
#include<stdlib.h>
#include"DDecayParameters.h"
#include"PhaseSpaceParameterisation.h"

int main(int argc, char *argv[]) {
  if(argc != 4) {
    return 0;
  }
  std::string filename_d = argv[1];
  std::string filename_cs = argv[2];
  std::string filename_K = argv[3];
  DDecayParameters ddecay(filename_d);
  ddecay.PlotParameters(filename_cs, filename_K);
  return 0;
}
