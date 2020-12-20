// Martin Duy Tat 30th October 2020
/**
 * EventList is a class that contains all events in a sample
 */

#ifndef EVENTLIST_H
#define EVENTLIST_H

#include<vector>
#include<string>
#include"Event.h"
#include"TTree.h"

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
    void AddEvent(const Event &event);
    /**
     * Function that returns total number of events in this EventList
     */
    int NumberEvents() const;
    /**
     * Function that returns the vector of Event objects
     * @return Vector of Event objects
     */
    const std::vector<Event>& GetEvents() const;
    /**
     * Function that returns an Event object
     * @param i Event index
     */
    const Event& GetEvent(const int &i) const;
    /**
     * Function for loading a TTree with events
     * @param filename Filename with TTree
     */
    void LoadTree(const std::string &filename);
    /**
     * Function for clearing the EventList to free up memory
     */
    void Clear();
  private:
    /**
     * List of events
     */
    std::vector<Event> m_eventlist;
};

#endif
