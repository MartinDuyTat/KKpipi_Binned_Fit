// Martin Duy Tat 12th December 2020
/**
 * KKpipiFit is a namespace with helper functions for doing the fitting and pull studies
 */

#ifndef KKPIPIFIT
#define KKPIPIFIT

#include<string>
#include"KKpipiFit.h"
#include"PhaseSpaceParameterisation.h"
#include"NaivePhaseSpace.h"
#include"RectangularPhaseSpace.h"
#include"SophisticatedPhaseSpace.h"
#include"AmplitudePhaseSpace.h"

namespace KKpipiFit {
  /**
   * Function for allocating the correct binning scheme
   * @param phasespace Pointer to a PhaseSpaceParameterisation object, where the chosen binning scheme is allocated
   */
  PhaseSpaceParameterisation* PickBinningScheme(const std::string &binning_choice);
}

#endif
