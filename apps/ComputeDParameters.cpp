// Martin Duy Tat 6th November 2020
/** 
 * This program is for calculating hadronic parameters of D->KKpipi decays using Monte Carlo integration and an amplitude model from AmpGen
 * @param 1 Choice of binning scheme, in lowercase letters
 * @param 2 Filename to save D hadronic decay parameters to
 * @param 3 Number of events used in Monte Carlo integration
 * @param 4 For Rectangular Phase Space, state the number of bins in each direction, for Sophisticated and Amplitude Phase Space state the total number of bins
 * @param 5 For Sophisticated Phase Space, input the filename for strong phases for binning scheme
 */

#include<string>
#include<stdlib.h>
#include<iostream>
#include"DDecayParameters.h"
#include"PhaseSpaceParameterisation.h"
#include"NaivePhaseSpace.h"
#include"RectangularPhaseSpace.h"
#include"SophisticatedPhaseSpace.h"
#include"AmplitudePhaseSpace.h"
#include"Constants.h"

int main(int argc, char *argv[]) {
  NaivePhaseSpace phasespace_naive;
  RectangularPhaseSpace phasespace_rectangular;
  SophisticatedPhaseSpace phasespace_sophisticated(atoi(argv[4]));
  AmplitudePhaseSpace phasespace_amplitude(atoi(argv[4]));
  PhaseSpaceParameterisation *psp;
  if(std::string(argv[1]) == "naive" && argc == 4) {
    psp = &phasespace_naive;
  } else if(std::string(argv[1]) == "rectangular" && argc == 9) {
    std::vector<int> bins = {atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8])};
    phasespace_rectangular = RectangularPhaseSpace(bins);
    psp = &phasespace_rectangular;
  } else if(std::string(argv[1]) == "sophisticated" && argc == 6) {
    psp = &phasespace_sophisticated;
    phasespace_sophisticated.ReadAverageStrongPhases(std::string(argv[5]));
  } else if(std::string(argv[1]) == "amplitude" && argc == 5) {
    psp = &phasespace_amplitude;
  } else {
    std::cout << "Invalid inputs!\n";
    return 0;
  }
  std::string filename = argv[2];
  int events = atoi(argv[3]);
  DDecayParameters ddecay(psp, events);
  ddecay.SaveCSV(filename);
  return 0;
}
