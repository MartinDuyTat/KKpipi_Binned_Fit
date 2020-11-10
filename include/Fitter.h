// Martin Duy Tat 1st November 2020
/**
 * Fitter is a class for maximising the likelihood and obtaining the CP violation parameters for a B meson decay
 */

#ifndef FITTER
#define FITTER

#include"BinList.h"
#include"DDecayParameters.h"
#include"CPParameters.h"

class Fitter {
  public:
    /**
     * Constructor that takes in a BinList object of input data and D meson decay parameters
     * @param binlist Input data events
     * @param ddecayparameters Parameters describing the D meson decay
     */
    Fitter(BinList binlist, DDecayParameters ddparameters);
    /**
     * Function for doing fit and returning the CP violation parameters (by reference)
     * @param cpparameters Initial guess of CP violation parameters, function replaces these with the fitted parameters and its covariance matrix
     */
    void DoFit(CPParameters &cpparameters);
  private:
    /**
     * Input data, sorted into their respective bins
     */
    BinList m_binlist;
    /**
     * D meson decay parameters
     */
    DDecayParameters m_ddparameters;
};

#endif
