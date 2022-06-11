// Martin Duy Tat 31st October 2020
/**
 * DDecayParameters is the class that calculates and stores the parameters describing the D^0 and DBAR^0 decay
 * These parameters only need to be calculated once because they only depend on the amplitude model
 */

#ifndef DDECAYPARAMETERS
#define DDECAYPARAMETERS

#include<vector>
#include<string>
#include"PhaseSpaceParameterisation.h"
#include"EventList.h"
#include"Amplitude.h"

class DDecayParameters {
  public:
    /**
     * Constructor that takes in a PhaseSpaceParameterisation object and calculates the D decay parameters in each bin
     * Unweighted events in phase space must be given in order to perform the Monte Carlo integrations
     * @param psp PhaseSpaceParameterisation object
     * @param eventlist EventList object with the unweighted Monte Carlo events
     * @param amplitude Optional argument, if given this Amplitude object will be used to calculate the strong phases
     */
    DDecayParameters(PhaseSpaceParameterisation *psp, const EventList &eventlist, const Amplitude *amplitude = nullptr);
    /** 
     * Constructor that takes in the D meson hadronic parameters from a comma separated CSV file, in the order i K_i Kbar_i c_i s_i
     * First line is assumed to be column names
     * @param filename Filename of file with D meson hadronic parameters
     */
    DDecayParameters(const std::string &filename);
    /**
     * Function for saving K_i, Kbar_i, c_i and s_i to a CSV file
     * @param filename Filename of file to save D meson hadronic parameters
     */
    void SaveCSV(const std::string &filename) const;
    /**
     * Function for plotting the c_i and s_i parameters in a s_i-c_i plane, and for plotting the magnitude of K_i
     * @param filename_cs Filename of c_i-s_i plot
     * @param filename_K Filename of K_i/Kbar_i plot
     */
    void PlotParameters(const std::string &filename_cs, const std::string &filename_K);
    /**
     * Function for getting fractional yield K_i
     * @return K Vector of fractional yields of D0 events
     */
    const std::vector<double>& GetK() const;
    /**
     * Function for getting fractional yield K_i
     * @return K Vector of fractional yields of DBAR0 events
     */
    const std::vector<double>& GetKbar() const;
    /**
     * Function for getting cosine of the strong phase
     * @return c Vector of cosine of the strong phases
     */
    const std::vector<double>& Getc() const;
    /**
     * Function for getting sine of the strong phase
     * @return s Vector of sine of the strong phases
     */
    const std::vector<double>& Gets() const;
    /**
     * Function for copying \f$c_i\f$ and \f$s_i\f$ from a DDecayParameter object with half the number of bins to this one
     * The \f$c_i\f$ and \f$s_i\f$ are repeated twice because the bins are constrained to have the same \f$c_i\f$ and \f$s_i\f$
     * The \f$K_i\f$ and \f$\bar{K_i}\f$ are kept the same
     */
    void CopyRedundantCS(const std::string &filename);
 private:
    /**
     * Vector of fractional yields of D0 events
     */
    std::vector<double> m_K;
    /**
     * Vector of fractional yields of DBAR0 events
     */
    std::vector<double> m_Kbar;
    /**
     * Vector of cosine of strong phases
     */
    std::vector<double> m_c;
    /**
     * Vector of sine of strong phases
     */
    std::vector<double> m_s;
    /**
     * Vector of bin areas of positive bins
     */
    std::vector<double> m_AreaPlus;
    /**
     * Vector of bin areas of negative bins
     */
    std::vector<double> m_AreaMinus;
};

#endif
