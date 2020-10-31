// Martin Duy Tat 31st October 2020
/**
 * BinList is a class that contains all the bins in phase space
 */

#ifndef BINLIST_H
#define BINLIST_H

#include<vector>
#include"PhaseSpaceParameterisation.h"
#include"Event.h"
#include"Bin.h"

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
  private:
    /**
     * Vector of Bin objects
     */
    std::vector<Bin> m_bins;
}/

#endif
