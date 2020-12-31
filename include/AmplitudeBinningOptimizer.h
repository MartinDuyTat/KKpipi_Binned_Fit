// Martin Duy Tat 31st December 2020
/**
 * AmplitudeOptimizer is a class that optimizes the AmplitudePhaseSpace binning by varying the bin edges along strong phase in order to maximize the binning \f$Q\f$-value
 */

#ifndef AMPLITUDEBINNINGOPTIMIZER
#define AMPLITUDEBINNINGOPTIMIZER

#include<string>
#include<vector>
#include"AmplitudePhaseSpace.h"
#include"EventList.h"

class AmplitudeBinningOptimizer {
  public:
    /**
     * Constructor that loads flat phase space events
     * @param bins Number of bins in the binning scheme, must be an even number
     * @param filename ROOT file containing flat phase space events
     * @param BinEdges Optional guess of initial bin edges, by default uniform bins are used for iniital guess
     */
  AmplitudeBinningOptimizer(const int &bins, const std::string filename, const std::vector<double> &BinEdges = std::vector<double>());
    /**
     * () operator overload that evaluates the \f$Q\f$-value of the binning scheme, given the bin edges
     * @param BinEdges Array of positive, increasing bin edges
     * @return Binning \f$-Q\f$-value
     */
    double operator ()(const double *BinEdges);
    /**
     * Function that runs the TMinuit2 fitting to maximize the \f$Q\f$-value
     * Iniital guess is a uniform binning, and optimized bin edges are saved in m_BinEdges
     */
    double OptimizeBinEdges();
    /**
     * Function that returns the bin edges
     * @return The optimized binning \f$Q\f$-value
     */
    const std::vector<double>& GetBinEdges() const;
  private:
    /**
     * Number of bins in the binning scheme
     */
    int m_bins;
    /**
     * Positive bin edges in increasing order
     * When the class is initialized, this variable contains the initial guess
     * After optimization the optimized bin edges are stored in this variable
     */
    std::vector<double> m_BinEdges;
    /**
     * EventList object with flat phase space events
     */
    EventList m_eventlist;
};

#endif
