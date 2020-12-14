// Martin Duy Tat 10th November 2020
/**
 * Gamma is a class that contains the parameters \f$r_B\f$, \f$\delta_B\f$ and \f$\gamma\f$ and their covariance matrix
 */

#ifndef GAMMA
#define GAMMA

#include"TMatrixD.h"

class Gamma {
  public:
    /**
     * Constructor that takes the CP parameters \f$r_B\f$, \f$\delta_B\f$ and \f$\gamma\f$
     * @param rB \f$r_B\f$
     * @param deltaB \f$\delta_B\f$
     * @param \f$\gamma\f$
     */
    Gamma(const double &rB, const double &deltaB, const double &gamma);
    /**
     * Function for getting CP parameters
     * @param rB \f$r_B\f$
     * @param deltaB \f$\delta_B\f$
     * @param \f$\gamma\f$
     */
    void GetGammaParameters(double &rB, double &deltaB, double &gamma) const;
    /**
     * Function for setting CP parameter covariance matrix
     * @param CovMatrix Array with covariance matrix
     */
    void SetCov(double *CovMatrix);
    /**
     * Function for getting CP parameter covariance matrix
     * @param CovMatrix TMatrixD object with covariance matrix
     */
    TMatrixD GetCov() const;
  private:
    /**
     * \f$r_B\r$
     */
    double m_rB;
    /**
     * \f$\delta_B\f$
     */
    double m_deltaB;
    /**
     * \f$\gamma\f$
     */
    double m_gamma;
    /**
     * Covariance matrix
     */
    TMatrixD m_CovMatrix;
};

#endif
