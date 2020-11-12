// Martin Duy Tat 12th November 2020
/**
 * RectangularPhaseSpace is a class that contains the information about how phase space is divided into bins using a rectangular parameterisation
 */

#ifndef RECTANGULARPHASESPACE
#define RECTANGULARPHASESPACE

#include"Event.h"
#include"PhaseSpaceParameterisation.h"

class RectangularPhaseSpace: public PhaseSpaceParameterisation {
  public:
    /**
     * Default constructor
     */
    RectangularPhaseSpace();
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
    int NumberOfBins() const;
  private:
};

#endif
