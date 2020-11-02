// Martin Duy Tat 2nd November 2020

#include<vector>
#include"EventList.h"
#include"Event.h"

EventList::Eventlist() {
}

void EventList::AddEvent(Event event) {
  m_eventlist.push_back(event);
}

int EventList::NumberEvents() {
  return m_eventlist.size();
}

std::vector<Event> EventList::GetEvents() {
  return m_eventlist;
}
