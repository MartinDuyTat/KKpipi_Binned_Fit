// Martin Duy Tat 30th October 2020
/**
 * Event is a class for storing the four-momenta of daugher particles in a D->KKpipi decay.
 */

#ifndef EVENT_H
#define EVENT_H

#include<vector>
#include"TLorentzVector.h"

class Event {
  public:
    /**
     * Default constructor for D to K+ K- pi+ pi- event with zero momentum
     */
    Event();
    /**
     * Constructor that takes a vector of four-momenta
     * @param p Four-momenta in the form (px, py, pz, E), in the order K+ K- pi+ pi-
     */
    Event(const std::vector<double> &p);
    /**
     * Returns the four-momenta of daughter particles as a vector
     * @return Four-momenta of daughter particles in the form (E, px, py, pz), in the order K+ K- pi+ pi-
     */
    const std::vector<double>& GetEventVector() const;
    /**
     * Constructor that takes a vector of four-momenta
     * @param p Vector of TLorentzVector objects, in the order K+ K- pi+ pi-
     */
    Event(const std::vector<TLorentzVector> &p);
    /**
     * Function for getting invariant mass of two particles
     * @param particle1 Particle 0(K+), 1(K-), 2(pi+), 3(pi-)
     * @param particle2 Particle 0(K+), 1(K-), 2(pi+), 3(pi-)
     * @return Returns invariant mass of given particles
     */
    double GetInvMass2(const int &particle1, const int &particle2) const;
    /**
     * Function for getting invariant mass of three particles
     * @param particle1 Particle 0(K+), 1(K-), 2(pi+), 3(pi-)
     * @param particle2 Particle 0(K+), 1(K-), 2(pi+), 3(pi-)
     * @param particle3 Particle 0(K+), 1(K-), 2(pi+), 3(pi-)
     * @return Returns mass of given particles
     */
    double GetInvMass3(const int &particle1, const int &particle2, const int &particle3) const;
  private:
    /**
     * Four-momenta of daughter particles
     */
    std::vector<double> m_momenta;
};

#endif
