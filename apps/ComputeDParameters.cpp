// Martin Duy Tat 6th November 2020
/** 
 * This program is for calculating hadronic parameters of D->KKpipi decays using Monte Carlo integration and an amplitude model from AmpGen
 * @param 1 Filename to save D hadronic decay parameters to
 * @param 2 Number of events used in Monte Carlo integration
 */

#include<string>
#include<stdlib.h>
#include"DDecayParameters.h"
#include"PhaseSpaceParameterisation.h"

int main(int argc, char *argv[]) {
  if(argc != 3) {
    return 0;
  }
  std::string filename = argv[1];
  int events = atoi(argv[2]);
  double mass_parent = 1.86483;
  double masses[4] = {0.493677, 0.493677, 0.13957039, 0.13957039};
  PhaseSpaceParameterisation psp;
  DDecayParameters ddecay(psp, mass_parent, masses, events);
  ddecay.SaveCSV(filename);
  return 0;
}
