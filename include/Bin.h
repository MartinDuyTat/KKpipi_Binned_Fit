// Martin Duy Tat 30th October 2020
/**
 * Bin is a class for a bin in phase space
 */
#ifndef BIN_H
#define BIN_H

#include"EventList.h"
#include"Event.h"

class Bin {
  public:
    /**
     * Default constructor that creates an empty EventList
     */
    Bin();
    /**
     * Function for adding an event
     * @param event Event to add
     * @param charge +1 for B+, -1 for B-
     */
    void AddEvent(Event event, int charge);
    /**
     * Function for getting number of events in this bin
     * @param charge +1 for B+, -1 for B-
     * @return Number of events in this bin
     */
    int GetNumberEvents(int charge);
  private:
    /**
     * EventList for B+ events in this bin
     */
    EventList m_eventlistBplus;
    /**
     * EventList for B- events in this bin
     */
    EventList m_eventlistBminus;
};

#endif
