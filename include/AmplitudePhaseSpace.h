// Martin Duy Tat 9th December 2020
/**
 * AmplitudePhaseSpace is a class that contains a binning scheme where the amplitude model is used to directly put events in bins
 */

#ifndef AMPLITUDEPHASESPACE
#define AMPLITUDEPHASESPACE

#include"PhaseSpaceParameterisation.h"
#include"Event.h"
#include"Amplitude.h"

class AmplitudePhaseSpace: public PhaseSpaceParameterisation {
  public:
    /**
     * Constructor that creates a binning with bins number of bins
     * @param bins Number of bins in the binning scheme
     */
    AmplitudePhaseSpace(const int &bins);
    /**
     * Virtual destructor to ensure well-defined behaviour when inheriting from PhaseSpaceParameterisation
     */
    ~AmplitudePhaseSpace();
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
    /**
     * Amplitude object to calculate the strong phase difference of an event
     */
    Amplitude m_amplitude;
};

#endif
