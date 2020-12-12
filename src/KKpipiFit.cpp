// Martin Duy Tat 12th December 2020

#include<iostream>
#include<string>
#include<vector>
#include"PhaseSpaceParameterisation.h"
#include"NaivePhaseSpace.h"
#include"RectangularPhaseSpace.h"
#include"SophisticatedPhaseSpace.h"
#include"AmplitudePhaseSpace.h"

namespace KKpipiFit {

  PhaseSpaceParameterisation* PickBinningScheme(const std::string &binning_choice) {
    PhaseSpaceParameterisation *phasespace;
    if(binning_choice == "naive") {
      std::cout << "Using Naive binning scheme\n";
      phasespace = new NaivePhaseSpace;
    } else if(binning_choice == "rectangular") {
      std::cout << "Using Rectangular binning scheme\n";
      std::vector<int> rectbins(5);
      std::cout << "Please input bins in each direction:\n";
      std::cin >> rectbins[0] >> rectbins[1] >> rectbins [2] >> rectbins[3] >> rectbins [4];
      phasespace = new RectangularPhaseSpace(rectbins);
    } else if(binning_choice == "sophisticated") {
      std::cout << "Using Sophisticated binning scheme\n";
      int bins;
      std::string phasefilename;
      std::cout << "Please input total number of bins:\n";
      std::cin >> bins;
      std::cout << "Please input file with mean strong phases:\n";
      std::cin >> phasefilename;
      phasespace = new SophisticatedPhaseSpace(bins, phasefilename);
    } else if(binning_choice == "amplitude") {
      std::cout << "Using Amplitude binning scheme\n";
      int bins;
      std::cout << "Please input total number of bins:\n";
      std::cin >> bins;
      phasespace = new AmplitudePhaseSpace(bins);
    } else {
      phasespace = nullptr;
    }
    return phasespace;
  }

}
