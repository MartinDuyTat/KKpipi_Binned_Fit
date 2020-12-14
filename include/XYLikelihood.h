// Martin Duy Tat 10th November 2020
/**
 * XYLikelihood is a class for calculating the likelihood of \f$x_\pm\f$, \f$y_\pm\f$, given the fitted values of \f$x_+\f$, \f$x_-\f$, \f$y_+\f$, \f$y_-\f$ and their covariance matrix, as a function of \f$r_B\f$, \f$\delta_B\f$ and \f$\gamma\f$
 * () operator is overloaded to make the likelihood function easily accessible
 */

#ifndef XYLIKELIHOOD
#define XYLIKELIHOOD

#include"TMatrixD.h"
#include"CPParameters.h"
#include"Gamma.h"

class XYLikelihood {
  public:
    /** 
     * Constructor that takes in the true values of \f$r_B\f$, \f$\delta_B\f$, \f$\gamma\f$, the fitted values of \f$x_+\f$, \f$x_-\f$, \f$y_+\f$, \f$y_-\f$ and their covariance matrix
     * @param cpparameters A CPParameters object containing the fitted \f$x_\pm\f$ and \f$y_\pm\f$ and the covariance matrix
     */
    XYLikelihood(const CPParameters &cpparameters);
    /**
     * Operator overload of () that returns -2 times the likelihood function
     * @param gamma An array of \f$r_B\f$, \f$\delta_B\f$ and \f$\gamma\f$, in degrees
     * @return -2ln(L), where L is the likelihood function
     */
    double operator()(const double *gamma);
  private:
    /**
     * Fitted value of \f$x_+\f$
     */
    double m_xplus;
    /**
     * Fitted value of \f$x_-\f$
     */
    double m_xminus;
    /**
     * Fitted value of \f$y_+\f$
     */
    double m_yplus;
    /**
     * Fitted value of \f$y_-\f$
     */
    double m_yminus;
    /** 
     * Inverse of the covariance matrix
     */
    TMatrixD m_InvCov;
};

#endif
