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
   * Mass of the \f$\K_S^0\f$ meson
   */
  const double MASS_KS = 0.497611;
  /**
   * \f$r_B\f$, as assumed in simulation
   */
  const double rB = 0.1;
  /**
   * \f$\delta_B\f$, as assumed in simulation, in radians
   */
  const double dB = 130.0*TMath::Pi()/180.0;
  /**
   * \f$\gamma\f$, as assumed in simulation, in radians
   */
  const double gamma = 75.0*TMath::Pi()/180.0;
  /**
   * \f$\delta_B\f$, as assumed in simulation, in degrees
   */
  const double dB_d = 130.0;
  /**
   * \f$\gamma\f$, as assumed in simulation, in degrees
   */
  const double gamma_d = 75.0;
  /**
   * \f$x_+\f$, as assumed in simulation
   */
  const double xplus = rB*TMath::Cos(dB + gamma);
  /**
   * \f$x_-\f$, as assumed in simulation
   */
  const double xminus = rB*TMath::Cos(dB - gamma);
  /**
   * \f$y_+\f$, as assumed in simulation
   */
  const double yplus = rB*TMath::Sin(dB + gamma);
  /**
   * \f$y_-\f$, as assumed in simulation
   */
  const double yminus = rB*TMath::Sin(dB - gamma);
}

#endif
