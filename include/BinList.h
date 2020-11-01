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
    BinList(const PhaseSpaceParameterisation &php);
    /**
     * Function for adding an event to the correct bin
     * @param event Event object to be added to the correct bin
     */
    void AddEvent(const Event &event);
    /**
     * Function for loading events from input data into their respective bins
     * @param tree A ROOT TTree in the AmpGen format containing all the input data events
     */
    void LoadTTree(const TTree &tree);
    /**
     * Function for getting the number of events in each bin
     * @return A vector of the number of events in each bin
     */
    std::vector<int> GetEvents();
    /**
     * Function for calculating the number of events in each bin, given the D decay parameters and the CP parameters
     * @param ddparameters A DDecayParameters object that describes the D meson decay
     * @param cpparameters A CPParameters object that describes the CP violation in the B meson decay
     * @param bins Number of bins
     * @return A vector of the predicted number of events in each bin
     */
    std::vector<int> Predict(const DDecayParameters &ddparameters, const CPParameters &cpparameters, const int &bins);
  private:
    /**
     * Vector of Bin objects
     */
    std::vector<Bin> m_bins;
}/

#endif
