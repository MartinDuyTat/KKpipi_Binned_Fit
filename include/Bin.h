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
     * Constructor that assigns a bin number
     * @param binnumber Bin number
     */
    Bin(int binnumber);
    /**
     * Constructor that assigns a bin number and an EventList
     * @param binnumber Bin number
     * @param eventlist Eventlist with events in this bin
     */
    Bin(int binnumber, const EventList &eventlist);
    /**
     * Function for adding an event
     * @param event Event to add
     */
    void AddEvent(const Event &event);
    /**
     * Function for getting number of events in this bin
     * @return Number of events in this bin
     */
    int GetNumberEvents();
  private:
    /**
     * EventList for events in this bin
     */
    EventList m_eventlist;
};

#endif
