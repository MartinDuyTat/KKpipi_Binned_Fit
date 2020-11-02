// Martin Duy Tat 30th October 2020
/**
 * EventList is a class that contains all events in a sample
 */

#ifndef EVENTLIST_H
#define EVENTLIST_H

#include"Event.h"
#include"TTree.h"
#include<vector>

class EventList {
  public:
    /**
     * Default constructor that creates an empty EventList
     */
    EventList();
    /**
     * Function that adds an Event to the EventList
     * @param event New Event object to be added to the EventList
     */
    void AddEvent(Event event);
    /**
     * Function that returns total number of events in this EventList
     */
    int NumberEvents();
  private:
    /**
     * List of events
     */
    std::vector<Event> m_eventlist;

#endif
