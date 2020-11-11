// Martin Duy Tat 10th November 2020
/**
 * FitGamma is a class for maximising the likelihood and obtaining the CP violation parameters \f$r_B\f$, \f$\delta_B\f$ and \f$\gamma\f$
 */

#ifndef FITGAMMA
#define FITGAMMA

#include<string>
#include"XYLikelihood.h"
#include"CPParameters.h"
#include"Minuit2/Minuit2Minimizer.h"
#include"Math/Functor.h"

class FitGamma {
  public:
    /**
     * Constructor that takes in a CPParameters object with fitted \f$x_\pm\f$ and \f$y_\pm\f$
     * @param cpparameters A CPParameters object with fitted \f$x_\pm\f$ and \f$y_\pm\f$ and their covariance matrix
     */
    FitGamma(CPParameters cpparameters);
    /**
     * Function for doing fit and returning the CP violation parameters (by reference)
     * @param cpparameters Initial guess of \f$r_B\f$, \f$\delta_B\f$ and \f$\gamma\f$
     */
    void DoFit(Gamma &GammaParams);
    /**
     * Function for plotting contours after fitting
     */
    void PlotContours(std::string Filename_rB_deltaB, std::string Filename_deltaB_gamma, std::string Filename_gamma_rB, unsigned int npoints) const;
    /**
     * Destructor to kill the minimiser
     */
    ~FitGamma();
  private:
    /**
     * Fitted \f$x_\pm\f$ and \f$y_\pm\f$
     */
    CPParameters m_cpparameters;
    /**
     * TMinuit2 Minimizer object
     */
    ROOT::Minuit2::Minuit2Minimizer *m_minimiser = nullptr;
    /**
     * Likelihood object 
     */
    XYLikelihood *m_likelihood = nullptr;
    /**
     * Likelihood functor
     */
    ROOT::Math::Functor m_fcn;
    
};

#endif
