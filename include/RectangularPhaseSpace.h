// Martin Duy Tat 12th November 2020
/**
 * RectangularPhaseSpace is a class that contains the information about how phase space is divided into bins using a rectangular parameterisation
 */

#ifndef RECTANGULARPHASESPACE
#define RECTANGULARPHASESPACE

#include<vector>
#include<map>
#include"Event.h"
#include"PhaseSpaceParameterisation.h"

class RectangularPhaseSpace: virtual public PhaseSpaceParameterisation {
  public:
    /**
     * Constructor that sets up the bins and phase space parameters
     * @param bins A vector with the number of bins in each direction, such that the total number of bins is the product
     * @param masses An array of the \f$(D^0, K^+, K^-, \pi^+, \pi^-)\f$ masses, with default the same as Ampgen's values
     */
    RectangularPhaseSpace(std::vector<int> bins);
    /**
     * Default constructor, puts the whole phase space into a single bin
     */
    RectangularPhaseSpace();
    /**
     * Virtual destructor to ensure well-defined behaviour when inheriting from PhaseSpaceParameterisation
     */
    ~RectangularPhaseSpace();
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
    /**
     * Function that returns the lower boundary of phase space
     * @param Coordinate to return the lower boundary for
     */
    double GetLowerBoundary(int coordinate) const;
    /**
     * Function that returns the upper boundary of phase space
     * @param Coordinate to return the upper boundary for
     */
    double GetUpperBoundary(int coordinate) const;
  private:
    /** 
     * The lower boundary of phase space
     */
    std::vector<double> m_xlow;
    /**
     * The upper boundary of phase space
     */
    std::vector<double> m_xhigh;
    /**
     * Vector with the number of bins in each direction
     */
    std::vector<int> m_Binning;
    /**
     * A vector of maps that connect the upper bin edges with the bin numbers
     */
    std::vector<std::map<double, int>> m_BinMap;
};

#endif
