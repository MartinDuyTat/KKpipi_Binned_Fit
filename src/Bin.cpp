// Martin Duy Tat 30th October 2020

#include"Bin.h"
#include"EventList.h"
#include"Event.h"

Bin::Bin() {
}

void Bin::AddEvent(Event event, int charge) {
  if(charge == + 1) {
    m_eventlistBplus.AddEvent(event);
  } else if(charge == -1) {
    m_eventlistBminus.AddEvent(event);
  }
}

int Bin::GetNumberEvents(int charge) {
  if(charge == +1) {
    return m_eventlistBplus.NumberEvents();
  } else if(charge == -1) {
    return m_eventlistBminus.NumberEvents();
  }
}

EventList Bin::GetEvents(int charge) {
  if(charge == +1) {
    return m_eventlistBplus;
  } else if(charge == -1) {
    return m_eventlistBminus;
  } else {
    return EventList();
  }
}
