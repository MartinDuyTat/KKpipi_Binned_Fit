// Martin Duy Tat 1st December 2020
/**
 * KKpipi_Constants is a namespace with the relevant constants in this analysis
 */

#ifndef KKPIPI_CONSTANTS
#define KKPIPI_CONSTANTS

#include"TMath.h"

namespace KKpipi_Constants {
  /**
   * Mass of the \f$D^0\f$ meson
   */
  const double MASS_D = 1.86483;
  /**
   * Mass of the \f$K^+\f$ meson
   */
  const double MASS_K = 0.493677;
  /**
   * Mass of the \f$\pi^+\f$ meson
   */
  const double MASS_PI = 0.13957039;
  /**
   * \f$r_B\f$, as assumed in simulation
   */
  const double rB = 0.1;
  /**
   * \f$\delta_B\f$, as assumed in simulation
   */
  const double dB = 130.0*TMath::Pi()/180.0;
  /**
   * \f$\gamma\f$, as assumed in simulation
   */
  const double gamma = 75.0*TMath::Pi()/180.0;
}

#endif
