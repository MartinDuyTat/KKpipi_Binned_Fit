// Martin Duy Tat 31st October 2020
/**
 * BinList is a class that contains all the bins in phase space
 * BinList also loads the input data and puts it in their respective bins
 */

#ifndef BINLIST_H
#define BINLIST_H

#include<vector>
#include"PhaseSpaceParameterisation.h"
#include"Event.h"
#include"Bin.h"
#include"DDecayParameters.h"
#include"CPParameters.h"

class Binlist {
  public;
    /**
     * Constructor that takes a PhaseSpaceParameterisation object and creates the bins
     * @param php A PhaseSpaceParameterisation object that defines the bins in the 5D phase space
     */
    BinList(PhaseSpaceParameterisation php);
    /**
     * Function for adding an event to the correct bin
     * @param event Event object to be added to the correct bin
     * @param charge +1 for B+, -1 for B-
     */
    void AddEvent(Event event, int charge);
    /**
     * Function for loading events from input data into their respective bins
     * @param tree A ROOT TTree in the AmpGen format containing all the input data events
     * @param charge +1 for B+, -1 for B-
     */
    void LoadTTree(const TTree &tree, int charge);
    /**
     * Function for getting the number of events in each bin
     * @param charge +1 for B+, -1 for B-
     * @return A vector of the number of events in each bin
     */
    std::vector<int> GetEvents(int charge);
    /**
     * Function for calculating the number of events in each bin, given the D decay parameters and the CP parameters
     * @param ddparameters A DDecayParameters object that describes the D meson decay
     * @param cpparameters A CPParameters object that describes the CP violation in the B meson decay
     * @param BplusEvents Vector of predicted number of B+ events
     * @param BminusEvents Vector of predicted number of B- events
     * @param totalBplus Total number of B+ events
     * @param totalBminus Total number of B- events
     */
    void Predict(const DDecayParameters &ddparameters, const CPParameters &cpparameters, std::vector<double> &BplusEvents, std::vector<double> &BminusEvents, int totalBplus, int totalBminus);
  private:
    /**
     * Vector of Bin objects
     */
    std::vector<Bin> m_bins;
    /**
     * A parameterisation of phase space
     */
    PhaseSpaceParameterisation m_php;
}/

#endif
