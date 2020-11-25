// Martin Duy Tat 25th November 2020
/**
 * SophisticatedPhaseSpace is a class that contains a binning scheme using the rectangular coordinates in RectangularPhaseSpace, but the binning is based on the strong phases obtained using the LHCb amplitude model
 */

#ifndef SOPHISTICATEDPHASESPACE
#define SOPHISTICATEDPHASESPACE

#include<vector>
#include<string>
#include"RectangularPhaseSpace.h"
#include"Event.h"
#include"Amplitude.h"

class SophisticatedPhaseSpace: public RectangularPhaseSpace {
  public:
    /**
     * Constructor that calls the corresponding constructor in RectangularPhaseSpace
     * @param bins A vector with the number of bins in each direction, such that the total number of bins is the product
     * @param masses An array of the \f$(D^0, K^+, K^-, \pi^+, \pi^-)\f$ masses, with default the same as Ampgen's values
     */
    SophisticatedPhaseSpace(std::vector<int> bins, double *masses = nullptr);
    /**
     * Default constructor that calls the corresponding constructor in RectangularPhaseSpace to put the whole phase space into a single bin
     */
    SophisticatedPhaseSpace();
    /**
     * Destructor that deletes the amplitude model
     */
    ~SophisticatedPhaseSpace();
    /**
     * Function for binning of the \f$(x_3, x_4)\f$ plane
     */
    void x3x4WhichBin(const Event &event, int &x3bin, int &x4bin) const;
    /**
     * Function for loading amplitude model
     */
    void LoadAmplitudeModel(const std::string &Damplitude, const std::string &DBARamplitude);
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
};

#endif
