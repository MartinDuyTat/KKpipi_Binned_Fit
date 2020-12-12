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
#include"CPParameters.h"
#include"Gamma.h"
#include"FitGamma.h"

namespace KKpipiFit {
  /**
   * Function for allocating the correct binning scheme
   * @param phasespace Pointer to a PhaseSpaceParameterisation object, where the chosen binning scheme is allocated
   */
  PhaseSpaceParameterisation* PickBinningScheme(const std::string &binning_choice);
  /**
   * Function for printing fit results for \f$x_\pm\f$ and \f$y_\pm\f$
   * @param cpparameters CPParameters object with fit results
   */
  void PrintXY (const CPParameters &cpparameters);
  /**
   * Function for printing fit results for \f$r_B\f$, \f$\delta_B\f$ and \f$\gamma\f$
   * @param gammaparams Gamma object with fit results
   */
  void PrintGamma (const Gamma &gammaparams);
  /**
   * Function that asks if user wants to draw contours of \f$r_B\f$, \f$\delta_B\f$ and \f$\gamma\f$ after fitting
   * If user inputs "yes", the contours are saved as PNG files
   * @param fitgamma FitGamma object with fit results
   */
  void DrawContours(const FitGamma &fitgamma);
}

#endif
