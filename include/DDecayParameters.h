// Martin Duy Tat 31st October 2020
/**
 * DDecayParameters is the class that calculates and stores the parameters describing the D^0 and DBAR^0 decay
 * These parameters only need to be calculated once because they only depend on the amplitude model
 */

#ifndef DDECAYPARAMETERS
#define DDECAYPARAMETERS

#include<vector>
#include"PhaseSpaceParameterisation.h"

class DDecayParameters {
  public:
    /**
     * Constructor that takes in a PhaseSpaceParameterisation object and calculates the D decay parameters in each bin
     * @param psp PhaseSpaceParameterisation object
     * @param events Number of events in each bin for Monte Carlo integration
     */
    DDecayParameters(const PhasespaceParameterisation &psp, const double &mass_parent, const double *mass_decay, int events);
    /**
     * Function for getting fractional yield K_i
     * @return K Vector of fractional yields of D0 events
     */
    std::vector<double> GetK();
    /**
     * Function for getting fractional yield K_i
     * @return K Vector of fractional yields of DBAR0 events
     */
    std::vector<double> GetKbar();
    /**
     * Function for getting cosine of the strong phase
     * @return c Vector of cosine of the strong phases
     */
    std::vector<double> Getc();
    /**
     * Function for getting sine of the strong phase
     * @return s Vector of sine of the strong phases
     */
    std::vector<double> Gets();
 private:
    /**
     * Vector of fractional yields of D0 events
     */
    std::vector<double> m_K;
    /**
     * Vector of fractional yields of DBAR0 events
     */
    std::vector<double> m_Kbar;
    /**
     * Vector of cosine of strong phases
     */
    std::vector<double> m_c;
    /**
     * Vector of sine of strong phases
     */
    std::vector<double> m_s;
}

#endif
