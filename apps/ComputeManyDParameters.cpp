// Martin Duy Tat 6th November 2020
/** 
 * This program is for calculating hadronic parameters of D->KKpipi decays using Monte Carlo integration and an amplitude model from AmpGen, and repeating the procedure many times
 * Unweighted events must be input in order to perform Monte Carlo integration
 * @param 1 Choice of binning scheme, in lowercase letters
 * @param 2 Filename prefix to save D hadronic decay parameters to
 * @param 3 Filename of TTree with unweighted events
 * @param 4 Filename prefix of \f$D^0\f$ model file
 * @param 5 Filename prefix of \f$\bar{D^0}\f$ model file
 * @param 6 Number of iterations
 * @param 7 First iteration number
 */

#include<string>
#include<stdlib.h>
#include<iostream>
#include"DDecayParameters.h"
#include"PhaseSpaceParameterisation.h"
#include"Constants.h"
#include"EventList.h"
#include"KKpipiFit.h"
#include"Amplitude.h"

int main(int argc, char *argv[]) {
  if(argc != 8) {
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
  eventlist.LoadTree(std::string(argv[3]));
  std::cout << "MC events ready\n";
  std::cout << "Calculating D hadronic parameters...\n";
  for(int i = std::atoi(argv[7]); i < std::atoi(argv[6]) + std::atoi(argv[7]); i++) {
    std::cout << "Iteration " << i << "\n";
    std::string D0Filename(argv[4]), D0barFilename(argv[5]);
    D0Filename += std::to_string(i) + ".so";
    D0barFilename += std::to_string(i) + ".so";
    Amplitude amplitude(D0Filename, D0barFilename);
    DDecayParameters ddecay(phasespace, eventlist, &amplitude);
    ddecay.SaveCSV(std::string(argv[2]) + std::to_string(i) + ".csv");
  }
  std::cout << "Calculation of D hadronic decay parameters complete\n";
  return 0;
}
