// Martin Duy Tat 19th December 2020
/**
 * GenerateFlat is the program for generating unweighted events in phase space using the hit-or-miss method
 * Event is saved to a TTree using the rectangular parameterisation
 * The \f$D^0\f$ and \f$\bar{D^0}\f$ amplitudes are also saved
 * D meson decay parameters are obtained from an input file
 * @param 1 Number of events to generate
 * @param 2 Filename of ROOT file to save TTree to
 */

#include<string>
#include<iostream>
#include<complex>
#include<stdlib.h>
#include"KKpipiMath.h"
#include"Constants.h"
#include"Generator.h"
#include"Amplitude.h"
#include"Event.h"
#include"TTree.h"
#include"TFile.h"

int main(int argc, char *argv[]) {
  if(argc != 3) {
    std::cout << "Incorrect number of inputs\n";
    return 0;
  }
  std::cout << "Starting generation of unweighted phase space events\n";
  double mass_decay[] = {KKpipi_Constants::MASS_K, KKpipi_Constants::MASS_K, KKpipi_Constants::MASS_PI, KKpipi_Constants::MASS_PI};
  Generator generator(KKpipi_Constants::MASS_D, mass_decay, 4);
  Amplitude amplitude;
  std::vector<double> X(5);
  std::complex<double> Damp;
  std::complex<double> DBARamp;
  TFile f(argv[2], "RECREATE");
  TTree tree("FlatEvents", "FlatEvents");
  tree.Branch("x1", X.data() + 0, "x1/D");
  tree.Branch("x2", X.data() + 1, "x2/D");
  tree.Branch("x3", X.data() + 2, "x3/D");
  tree.Branch("x4", X.data() + 3, "x4/D");
  tree.Branch("x5", X.data() + 4, "x5/D");
  tree.Branch("D_amplitude", &Damp);
  tree.Branch("Dbar_amplitude", &DBARamp);
  std::cout << "TTree initialized\n";
  for(int i = 0; i < atoi(argv[1]); i++) {
    if(i%100000 == 1) {
      std::cout << i - 1 << " events generated\n";
    }
    Event event(generator.Generate());
    std::vector<double> coordinates = KKpipiMath::RectCoordinates(event);
    X = coordinates;
    Damp = amplitude(event.GetEventVector(), +1);
    DBARamp = amplitude(event.GetEventVector(), -1);
    tree.Fill();
  }
  std::cout << "Event generation finished, generated " << atoi(argv[1]) << " events in total\n";
  tree.Write();
  f.Close();
  std::cout << "Events saved to file\n";
  return 0;
}
