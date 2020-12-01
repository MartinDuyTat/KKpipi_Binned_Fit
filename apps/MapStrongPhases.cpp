// Martin Duy Tat 1st December 2020
/**
 * MapStrongPhases is a program that will calculate strong phases on a grid in phase space and save it to a file
 */

#include<iostream>
#include"KKpipiMath.h"
#include"TRandom3.h"

int main(int argc, char *argv[]) {
  if(argc != 6) {
    std::cout << "Need 5 inputs\n";
    return 0;
  }
  std::vector<double> X = {atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5])};
  std::vector<double> four_momenta = KKpipiMath::ConvertXToMomenta(X);
  std::vector<double> event = KKpipiMath::RectCoordinates(Event(four_momenta));
  std::cout << "Input coordinates:      " << X[0] << " " << X[1] << " " << X[2] << " " << X[3] << " " << X[4] << std::endl;
  std::cout << "Calculated coordinates: " << event[0] << " " << event[1] << " " << event [2] << " " << event[3] << " " << event[4] << std::endl;
  std::cout << four_momenta[0] + four_momenta[4] + four_momenta[8] + four_momenta[12] << std::endl;
  std::cout << four_momenta[1] + four_momenta[5] + four_momenta[9] + four_momenta[13] << std::endl;
  std::cout << four_momenta[2] + four_momenta[6] + four_momenta[10] + four_momenta[14] << std::endl;
  std::cout << four_momenta[3] + four_momenta[7] + four_momenta[11] + four_momenta[15] << std::endl;
  return 0;
}
