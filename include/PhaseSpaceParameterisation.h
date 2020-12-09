// Martin Duy Tat 12th November 2020
/**
 * PhaseSpaceParameterisation is a class that contains the information about how phase space is divided into bins
 * Add a new binning scheme by adding a derived class with the binning implemented
 */

#ifndef PHASESPACEPARAMETERISATION
#define PHASESPACEPARAMETERISATION

#include"Event.h"

class PhaseSpaceParameterisation {
  public:
    /**
     * Constructor that creates a binning scheme with bins number of bins
     * @param bins Number of bins in the binning scheme
     */
    PhaseSpaceParameterisation(int bins);
    /**
     * Function that determines which bin an event belongs to
     * @param event The event we want to determine the bin of
     * @return Bin number
     */
    virtual int WhichBin(const Event &event) const = 0;
    /**
     * Function that returns the number of bins in the binning scheme
     * @return Number of bins
     */
    int NumberOfBins() const;
  private:
    /**
     * Number of bins in this binning scheme
     */
    int m_bins;
};

#endif
