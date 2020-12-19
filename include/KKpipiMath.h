// Martin Duy Tat 1st December 2020
/**
 * KKpipiMath is a namespace that contain all the maths functions used in the KKpipi library
 */

#ifndef KKPIPIMATH
#define KKPIPIMATH

#include<vector>
#include"TLorentzVector.h"
#include"Event.h"
#include"CPParameters.h"
#include"DDecayParameters.h"

namespace KKpipiMath {
  /**
   * Function that converts a vector into a vector of TLorentzVector objects
   * @param event A vector with four-momenta, in the order \f$(K^+, K^-, \pi^+, \pi^-)\f$
   * @return A vector of TLorentzVector objects
   */
  std::vector<TLorentzVector> ConvertEventTo4Vectors(const std::vector<double> &event);
  /**
   * Function that determines the rectangular coordinates of an event
   * @param event Vector of four-momenta of the final states in the event, in the order \f$(K^+, K^-, \pi^+, \pi^-)\f$
   * @return A vector with the five rectangular coordinates
   */
  std::vector<double> RectCoordinates(const Event &event);
  /**
   * Function that converts a set of five rectangular coordinates into a vector of four-momenta
   * @param X Vector with five rectangular coordinates
   * @return Vector with four-momenta of the event, in the order \f$(K^+, K^-, \pi^+, \pi^-)\f$
   */
  std::vector<double> ConvertXToMomenta(const std::vector<double> &X);
  /**
   * Function that predicts how many events there are in each bin, given \f$K_i\f$, \f$\bar{K_i}\f$, \f$c_i\f$, \f$s_i\f$, \f$x_\pm\f$, \f$y_\pm\f$ and the total number of \f$B^\pm$ events
   * @param ddparameters A DDecayParameters object that describes the D meson decay
   * @param cpparameters A CPParameters object that describes the CP violation in the B meson decay
   * @param BplusEvents Vector of predicted number of B+ events
   * @param BminusEvents Vector of predicted number of B- events
   * @param totalBplus Total number of B+ events
   * @param totalBminus Total number of B- events
   */
  void ExpectedNumberOfEvents(const DDecayParameters &ddparameters,
	       const CPParameters &cpparameters,
	       const int &totalBplus,
	       const int &totalBminus,
	       std::vector<double> &BplusEvents,
	       std::vector<double> &BminusEvents);
}

#endif
