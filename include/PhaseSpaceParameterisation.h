// Martin Duy Tat 12th November 2020
/**
 * PhaseSpaceParameterisation is a class that contains the information about how phase space is divided into bins
 * Add a new binning scheme by adding a derived class with the binning implemented
 * Bin numbers go from \f$i = 1, 2, 3, ...\f$, and \f$i = -1, -2, -3, ...\f$ for CP conjugated bins
 * Which events that are CP conjugates are defined by the binning scheme
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
    PhaseSpaceParameterisation(const int &bins);
    /**
     * Virtual destructor for well-defined behaviour when using polymorphism
     */
    virtual ~PhaseSpaceParameterisation() = 0;
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
    /**
     * Function that sets the veto boundary
     */
    void SetKSVeto(double veto);
    /**
     * Function that returns true if event is inside the veto window
     */
    bool isKSVeto(const Event &event) const;
  private:
    /**
     * Number of bins in this binning scheme
     */
    int m_bins;
    /**
     * Veto window size around the \f$K_S^0\f$ mass, in units of \f$\text{GeV}\f$
     * If this is negative or zero, no veto is applied
     */
    double m_KSVeto = 0.0;
};

#endif
