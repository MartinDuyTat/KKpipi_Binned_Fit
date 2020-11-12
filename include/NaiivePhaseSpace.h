// Martin Duy Tat 31st October 2020
/**
 * NaiivePhaseSpace is a class that contains the information about how phase space is divided into bins in a naiive way
 * NaiivePhaseSpace contains a very coarse and arbitrary binning of phase space
 */

#ifndef NAIIVEPHASESPACE
#define NAIIVEPHASESPACE

#include"Event.h"
#include"PhaseSpaceParameterisation.h"

class NaiivePhaseSpace: public PhaseSpaceParameterisation {
  public:
    /**
     * Default constructor that sets the number of bins
     */
    NaiivePhaseSpace();
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
