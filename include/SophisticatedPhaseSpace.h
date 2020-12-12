// Martin Duy Tat 25th November 2020
/**
 * SophisticatedPhaseSpace is a class that contains a binning scheme using the rectangular coordinates in RectangularPhaseSpace, but the binning is based on the strong phases obtained using the LHCb amplitude model
 */

#ifndef SOPHISTICATEDPHASESPACE
#define SOPHISTICATEDPHASESPACE

#include<vector>
#include<string>
#include<map>
#include"RectangularPhaseSpace.h"
#include"PhaseSpaceParameterisation.h"
#include"Event.h"
#include"Amplitude.h"

class SophisticatedPhaseSpace: virtual public PhaseSpaceParameterisation, public RectangularPhaseSpace {
  public:
    /**
     * Constructor that calls the corresponding constructor in RectangularPhaseSpace
     * @param nbins Number of bins in total
     */
    SophisticatedPhaseSpace(int nbins);
    /**
     * Function for reading in average strong phases in the \f$(x_1, x_2, x_5)\f$ volume from a file and storing them in a vector
     * @param filename Filename of file with mean strong phases
     */
    void ReadAverageStrongPhases(const std::string &filename);
    /**
     * Constructor that also reads in the average strong phases and sets up a lookup table
     */
    SophisticatedPhaseSpace(int nbins, const std::string &filename);
    /**
     * Destructor that deletes the amplitude model
     * Virtual destructor to ensure well-defined behaviour when inheriting from PhaseSpaceParameterisation
     */
    ~SophisticatedPhaseSpace();
    /**
     * Function for binning of the \f$(x_3, x_4)\f$ plane on a grid
     */
    void x3x4WhichBinGrid(const Event &event, int &x3bin, int &x4bin) const;
    /**
     * Function for loading amplitude model
     */
    void LoadAmplitudeModel();
    /**
     * Function for calculating the strong phase difference of an event
     */
    double StrongPhase(const Event &event) const;
    /**
     * Function for calculating the average strong phase and its standard deviation in each bin in the \f$(x_3, x_4)\f$ plane
     * @param BplusFilename Filename with ROOT file with \f$B^+\f$ events in AmpGen format
     * @param BminusFilename Filename with ROOT file with \f$B^-\f$ events in AmpGen format
     * @param MeanFilename Filename of CSV file that contains the mean strong phases
     * @param RMSFilename Filename of CSV file that contains the RMS of phases
     */
    void CalculateStrongPhases(std::string BplusFilename, std::string BminusFilename, std::string MeanFilename, std::string RMSFilename) const;
    /**
     * Function for clearing out vector of phases to free up memory
     */
    void ClearAverageStrongPhases();
    /**
     * Function that returns the number of regions in the \f$(x_3, x_4)\f$ plane
     */
    int NumberOfRegions() const;
    /**
     * Function that determines which region in the \f$(x_3, x_4)\f$ plane an event belongs to
     * @param event Vector with four-momenta of the event
     * @return Region number
     */
    int WhichRegion(const std::vector<double> &X) const;
    /**
     * Function that determines which bin an event belongs to, based on the lookup table
     * The first m_binregion regions are separate bins and assinged to the first m_binregion bins without using a lookup table
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
     * Number of regions in the \f$(x_3, x_4)\f$ plane
     */
    int m_regions;
    /**
     * Number of regions that are also separate bins
     */
    int m_binregion;
    /**
     * Number of bins in the \f$x_3\f$ direction
     */
    int m_x3bins = 100;
    /**
     * Number of bins in the \f$x_4\f$ direction
     */
    int m_x4bins = 100;
    /**
     * Amplitude object to calculate event amplitudes
     */
    Amplitude *m_AmplitudeModel = nullptr;
    /**
     * 4-dimensional lookup vector with bin numbers
     * Indices are: bin number x_1 x_2 x_5
     */
    std::vector<std::vector<std::vector<std::vector<int>>>> m_LookupBins;
};

#endif
