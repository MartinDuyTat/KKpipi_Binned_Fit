// Martin Duy Tat 9th December 2020
/**
 * AmplitudePhaseSpace is a class that contains a binning scheme where the amplitude model is used to directly put events in bins
 */

#ifndef AMPLITUDEPHASESPACE
#define AMPLITUDEPHASESPACE

#include<map>
#include<vector>
#include"PhaseSpaceParameterisation.h"
#include"Event.h"
#include"Amplitude.h"

class AmplitudePhaseSpace: public PhaseSpaceParameterisation {
  public:
    /**
     * Constructor that creates a binning with bins number of bins
     * If using bins along \f$r_D\f$, the total number of bins must be even!
     * @param bins Number of bins in the binning scheme
     * @param rDBinning True when using \f$4\f$ bins in the \f$r_D\f$ direction, instead of simply splitting at \f$r_D = 1\f$
     * @param D0Filename Load alternative model by providing a filename
     * @param Dbar0Filename Load alternative model by providing a filename
     */
    AmplitudePhaseSpace(const int &bins, bool rDBinning = false, const std::string &D0Filename = "", const std::string &Dbar0Filename = "");
    /**
     * Virtual destructor to ensure well-defined behaviour when inheriting from PhaseSpaceParameterisation
     */
    ~AmplitudePhaseSpace();
    /**
     * Function for setting the flag for reading amplitudes from the Event object
     * Set true if amplitude is read from the Event Object
     * @param TrueIfReadFromEvent Obvious name
     */
    void ReadAmplitudeFromEvent(bool TrueIfReadFromEvent);
    /**
     * Function for setting the flag m_UseVariableBinWidths
     * Before setting this flag to true, the map BinEdgeMap must be initialised with the correct number of bin edges, otherwise this flag stays false
     * @param VariableBinWidth Flag that is true when variable bin widths are used
     */
    void UseVariableBinWidths(bool VariableBinWidths);
    /**
     * Function for setting the bin edges when variable bin widths are used
     * For \f$N\f$ bins, \f$(N - 2)/2\f$ inputs are required
     * \f$\delta_D = 0, \pm\pi\f$ are fixed bin edges, and the other $\f$N - 2\f$ are placed symmetrically around $\f$\delta = 0\f$, therefore $\f$(N - 2)/2\f$ bin edges must be specified
     * Bin edges must be positive, in increasing order!
     */
    void SetBinEdges(const std::vector<double> &BinEdges);
    /**
     * Function that determines which bin an event belongs to
     * @param event The event we want to determine the bin of
     * @return Bin number
     */
    int WhichBin(const Event &event) const;
    /**
     * Function that returns the number of bins in the binning scheme
     * When using \f$4\f$ bins along the \f$r_D\f$ direction, this returns the number of bins divided by two because practically half of the bins are redundant (however, this function is not virtual so calling it from a parent pointer will return the total number of bins)
     * @return Number of bins
     */
    int NumberOfBins() const;
  private:
    /**
     * Amplitude object to calculate the strong phase difference of an event
     */
    Amplitude m_amplitude;
    /**
     * When true, the amplitudes are read from the Event object instead of calculated
     * This only works for MC events because the amplitude has already been calculated and stored
     */
    bool m_ReadAmplitudeFromEvent = false;
    /**
     * When true, the bin edges will have variable bin widths along strong phase difference instead of uniformly spaced
     */
    bool m_UseVariableBinWidths = false;
    /**
     * Map that connects the upper bin edges with the bin numbers
     */
    std::map<double, int> m_BinMap;
    /**
     * Flag for binning along \f$r_D\f$, in addition to binning in strong phase
     * If true, \f$4\f$ bins will be used, with bin boundaries at $\f$\ln(r_D) = 0\f$ and \f$\ln(r_D) = \pm 0.5\f$
     * If false, \f$2\f" bins will be used, split along $\f$r_D = 1\f$
     */
    bool m_rDBinning;
};

#endif
