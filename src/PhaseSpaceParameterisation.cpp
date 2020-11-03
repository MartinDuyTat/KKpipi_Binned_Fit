// Martin Duy Tat 2nd November 2020

#include"PhaseSpaceParameterisation.h"
#include"Event.h"

PhaseSpaceParameterisation::PhaseSpaceParameterisation(): m_bins(8) {
}

int PhaseSpaceParameterisation::WhichBin(const Event &event) {
  double s02 = event.GetInvMass2(0, 2);
  double s13 = event.GetInvMass2(1, 3);
  double s012 = event.GetInvMass3(0, 1, 2);
  if(s02 < 0.8 && s13 < 0.8 && s012 < 1.5) {
    return 0;
  } else if(s02 >= 0.8 && s13 < 0.8 && s012 < 1.5) {
    return 1;
  } else if(s02 < 0.8 && s13 >= 0.8 && s012 < 1.5) {
    return 2; 
  } else if(s02 < 0.8 && s13 < 0.8 && s012 >= 1.5) {
    return 3;
  } else if(s02 >= 0.8 && s13 >= 0.8 && s012 < 1.5) {
    return 4;
  } else if(s02 >= 0.8 && s13 < 0.8 && s012 >= 1.5) {
    return 5;
  } else if(s02 < 0.8 && s13 >= 0.8 && s012 >= 1.5) {
    return 6;
  } else {
    return 7;
  }
}

int PhaseSpaceParameterisation::NumberOfBins() {
  return m_bins;
}
