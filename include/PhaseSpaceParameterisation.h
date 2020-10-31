// Martin Duy Tat 31st October 2020
/**
 * PhaseSpaceParameterisation is a class that contains the information about how phase space is divided into bins
 * PhasespaceParameterisation contains a very coarse and arbitrary binning of phase space
 * A more sophisticated binning can be added by added a new class that inherits from PhaseSpaceParameterisation
 */

#ifndef PHASESPACEPARAMETERISATION
#define PHASESPACEPARAMETERISATION

#include"Event.h"

class PhaseSpaceParameterisation {
  public:
    /**
     * Default constructor
     */
    PhaseSpaceParameterisation();
    /**
     * Function that determines which bin an event belongs to
     * @param event The event we want to determine the bin of
     * @return Bin number
     */
    int WhichBin(const Event &event);
    /**
     * Function that returns the number of bins in the binning scheme
     * @return Number of bins
     */
    int NumberOfBins();
  private:
    /**
     * Number of bins in this binning scheme
     */
    int m_bins;
