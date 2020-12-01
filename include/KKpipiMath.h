// Martin Duy Tat 1st December 2020
/**
 * KKpipiMath is a namespace that contain all the maths functions used in the KKpipi library
 */

#ifndef KKPIPIMATH
#define KKPIPIMATH

#include<vector>
#include"TLorentzVector.h"
#include"Event.h"

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
}

#endif
