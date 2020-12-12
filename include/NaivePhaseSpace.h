// Martin Duy Tat 31st October 2020
/**
 * NaivePhaseSpace is a class that contains the information about how phase space is divided into bins in a naiive way
 * NaivePhaseSpace contains a very coarse and arbitrary binning of phase space
 */

#ifndef NAIVEPHASESPACE
#define NAIVEPHASESPACE

#include"Event.h"
#include"PhaseSpaceParameterisation.h"

class NaivePhaseSpace: public PhaseSpaceParameterisation {
  public:
    /**
     * Default constructor that sets the number of bins
     */
    NaivePhaseSpace();
    /**
     * Virtual destructor to ensure well-defined behaviour when inheriting from PhaseSpaceParameterisation
     */
    ~NaivePhaseSpace();
    /**
     * Function that determines which bin an event belongs to
     * @param event The event we want to determine the bin of
     * @return Bin number
     */
    int WhichBin(const Event &event) const;
    /**
     * Function that returns the number of bins in the binning scheme
     * @return Number of bins
     */
    int NumberOfBins() const;
  private:
};

#endif
