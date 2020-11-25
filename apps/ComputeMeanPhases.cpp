// Martin Duy Tat 25th November 2020

#include<iostream>
#include<string>
#include"SophisticatedPhaseSpace.h"

int main(int argc, char *argv[]) {
  if(argc != 5) {
    std::cout << "Incorrect number of inputs!\n";
    return 0;
  }
  SophisticatedPhaseSpace phasespace;
  phasespace.LoadAmplitudeModel("D0toKKpipi.so", "Dbar0toKKpipi.so");
  std::string BplusFilename = argv[1];
  std::string BminusFilename = argv[2];
  std::string MeanFilename = argv[3];
  std::string RMSFilename = argv[4];
  phasespace.CalculateStrongPhases(BplusFilename, BminusFilename, MeanFilename, RMSFilename);
  return 0;
}
