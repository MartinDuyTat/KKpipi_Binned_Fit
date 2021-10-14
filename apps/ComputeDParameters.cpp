// Martin Duy Tat 6th November 2020
/** 
 * This program is for calculating hadronic parameters of D->KKpipi decays using Monte Carlo integration and an amplitude model from AmpGen
 * Unweighted events must be input in order to perform Monte Carlo integration
 * @param 1 Choice of binning scheme, in lowercase letters
 * @param 2 Filename to save D hadronic decay parameters to
 * @param 3 Filename of TTree with unweighted events
 * @param 4 Filename of TTree with efficiency weights (optional)
 */

#include<string>
#include<stdlib.h>
#include<iostream>
#include"DDecayParameters.h"
#include"PhaseSpaceParameterisation.h"
#include"Constants.h"
#include"EventList.h"
#include"KKpipiFit.h"

int main(int argc, char *argv[]) {
  if(argc != 4 && argc != 5) {
    std::cout << "Incorrect number of inputs\n";
  }
  std::cout << "Calculation of K_i, Kbar_i, c_i, s_i using Monte Carlo integration\n";
  std::cout << "Loading phase space binning...\n";
  PhaseSpaceParameterisation *phasespace = KKpipiFit::PickBinningScheme(std::string(argv[1]));
  if(phasespace == nullptr) {
    std::cout << "Invalid binning scheme\n";
    return 0;
  }
  std::cout << "Binning scheme prepared\n";
  std::cout << "Loading MC events...\n";
  EventList eventlist;
  argc == 4 ? eventlist.LoadTree(std::string(argv[3])) : eventlist.LoadTree(std::string(argv[3]), std::string(argv[4]));
  std::cout << "MC events ready\n";
  std::cout << "Calculating D hadronic parameters...\n";
  DDecayParameters ddecay(phasespace, eventlist);
  ddecay.SaveCSV(std::string(argv[2]));
  std::cout << "Calculation of D hadronic decay parameters complete\n";
  return 0;
}
