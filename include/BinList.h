// Martin Duy Tat 31st October 2020
/**
 * BinList is a class that contains all the bins in phase space
 * BinList also loads the input data and puts it in their respective bins
 * Events can also be saved, instead of simply counted, by setting SaveEvents = true
 * Bin \f$-i\f$ is the CP conjugation of bin \f$i\f$, defined by the binning scheme through the function PhaseSpaceParameterisation::WhichBin
 */

#ifndef BINLIST_H
#define BINLIST_H

#include<vector>
#include"PhaseSpaceParameterisation.h"
#include"Event.h"
#include"Bin.h"
#include"DDecayParameters.h"
#include"CPParameters.h"

class BinList {
  public:
    /**
     * Constructor that takes a PhaseSpaceParameterisation object and creates the bins
     * @param php A PhaseSpaceParameterisation object that defines the bins in the 5D phase space
     */
    BinList(PhaseSpaceParameterisation *php, bool SaveEvents = false);
    /**
     * Function for adding an event to the correct bin
     * @param event Event object to be added to the correct bin
     * @param charge +1 for B+, -1 for B-
     */
    void AddEvent(const Event &event, const int &charge);
    /**
     * Function for adding an event to the correct bin, if the number of events in that bin is less than the maximum
     * @param event Event object to be added to the correct bin
     * @param charge +1 for B+, -1 for B-
     * @param maxEvents Maximum number of events in each bin
     */
    void AddEvent(const Event &event, const int &charge, const int &maxevents);
    /**
     * Function for loading events from input data into their respective bins
     * @param tree A ROOT TTree in the AmpGen format containing all the input data events
     * @param charge +1 for B+, -1 for B-
     * @param StartEvent Index of first event to read from TTree
     * @param TotalEvents Total number of events to read from file, -1 means read until end of file
     */
    void LoadTTree(TTree *tree, const int &charge, const int &StartEvent = 0, const int &TotalEvents = -1);
    /**
     * Function for getting the total number of bins, not counting CP conjugated bins
     */
    int NumberBins();
    /**
     * Function for getting the number of events in each bin
     * @param charge +1 for B+, -1 for B-
     * @param cp \f$+1\f$ for positive bins, \f$-1\f$ for negative bins
     * @return A vector of the number of events in each bin
     */
    std::vector<int> GetEvents(const int &charge, const int &cp) const;
    /**
     * Function for getting Bin object
     * @param i Bin number
     * @param cp \f$+1\f$ for positive bins, \f$-1\f$ for negative bins
     * @return Bin object
     */
    Bin GetBin(const int &i, const int &cp);
  private:
    /**
     * A parameterisation of phase space
     */
    PhaseSpaceParameterisation *m_psp;
    /**
     * Vector of Bin objects
     */
    std::vector<Bin> m_bins;
    /**
     * Vector of CP conjugated Bin objects
     */
    std::vector<Bin> m_CPbins;
    /**
     * Flag to save all event data into bins for later analysis, instead of simply counting the number of events
     */
    bool m_SaveEvents;
    /**
     * Vector containing the number of \f$B^+\f$ events in each bin
     */
    std::vector<int> m_BplusEvents;
    /**
     * Vector containing the number of \f$B^-\f$ events in each bin
     */
    std::vector<int> m_BminusEvents;
    /**
     * Vector containing the number of CP conjugated \f$B^+\f$ events in each bin
     */
    std::vector<int> m_BplusCPEvents;
    /**
     * Vector containing the number of CP conjugated \f$B^-\f$ events in each bin
     */
    std::vector<int> m_BminusCPEvents;
};

#endif
