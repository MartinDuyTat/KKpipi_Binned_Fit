// Martin Duy Tat 30th October 2020

#include"Bin.h"
#include"EventList.h"
#include"Event.h"

Bin::Bin() {
}

void Bin::AddEvent(const Event &event, const int &charge) {
  if(charge == + 1) {
    m_eventlistBplus.AddEvent(event);
  } else if(charge == -1) {
    m_eventlistBminus.AddEvent(event);
  }
}

int Bin::GetNumberEvents(const int &charge) const {
  if(charge == +1) {
    return m_eventlistBplus.NumberEvents();
  } else if(charge == -1) {
    return m_eventlistBminus.NumberEvents();
  }
  return 0;
}

EventList Bin::GetEvents(const int &charge) {
  if(charge == +1) {
    return m_eventlistBplus;
  } else if(charge == -1) {
    return m_eventlistBminus;
  } else {
    return EventList();
  }
}
