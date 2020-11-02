// Martin Duy Tat 31st October 2020
/**
 * Likelihood is a class for calculating the likelihood, given an EventList of input data and a set of D meson decay parameters and CP violation parameters in B meson decays
 * () operator is overloaded to make the likelihood function easily accessible
 */

#ifndef LIKELIHOOD
#define LIKELIHOOD

#include"BinList.h"
#include"DDecayParameters.h"

class Likelihood {
  public:
    /**
     * Constructor that takes in an BinList object with input data and a DDecayParameters object
     * @param events BinList object with the input data
     * @param ddecayparameters A DDecayParameters object with the parameters for the D meson decay
     */
    Likelihood(BinList bins, DDecayParameters ddparameters);
    /**
     * Operator overload of () to easily access the likelihood function
     * @param cpparameters A CPParameters object with the CP violation parameters for the B meson decay
     * @return -2*ln(L), where L is the likelihood function
     */
    double operator()(const double *cpparameters);
  private:
    /**
     * List of input data events
     */
    BinList m_bins;
    /**
     * D meson decay parameters
     */
    DDecayParameters m_ddparameters;
};

#endif
