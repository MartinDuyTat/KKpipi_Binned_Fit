// Martin Duy Tat 6th November 2020
/** 
 * This program is for calculating hadronic parameters of D->KKpipi decays using Monte Carlo integration and an amplitude model from AmpGen
 * Unweighted events must be input in order to perform Monte Carlo integration
 * @param 1 Choice of binning scheme, in lowercase letters
 * @param 2 Filename to save D hadronic decay parameters to
 * @param 3 Filename of TTree with unweighted events
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
#include"EventList.h"

int main(int argc, char *argv[]) {
  NaivePhaseSpace phasespace_naive;
  RectangularPhaseSpace phasespace_rectangular;
  SophisticatedPhaseSpace phasespace_sophisticated(atoi(argv[4]));
  AmplitudePhaseSpace phasespace_amplitude(atoi(argv[4]));
  PhaseSpaceParameterisation *psp;
  if(std::string(argv[1]) == "naive" && argc == 5) { // TODO Use choice function for binning scheme
    psp = &phasespace_naive;
  } else if(std::string(argv[1]) == "rectangular" && argc == 9) {
    std::vector<int> bins = {atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), atoi(argv[7]), atoi(argv[8])};
    phasespace_rectangular = RectangularPhaseSpace(bins);
    psp = &phasespace_rectangular;
  } else if(std::string(argv[1]) == "sophisticated" && argc == 6) {
    psp = &phasespace_sophisticated;
    phasespace_sophisticated.ReadAverageStrongPhases(std::string(argv[5]));
  } else if(std::string(argv[1]) == "amplitude" && argc == 5) {
    phasespace_amplitude.ReadAmplitudeFromEvent(true);
    psp = &phasespace_amplitude;
  } else {
    std::cout << "Invalid inputs!\n";
    return 0;
  }
  std::string filename = argv[2];
  EventList eventlist;
  eventlist.LoadTree(std::string(argv[3]));
  DDecayParameters ddecay(psp, eventlist);
  ddecay.SaveCSV(filename);
  return 0;
}
