// Martin Duy Tat 31st October 2020
/**
 * Generator is a class that generates uniformly distributed events in phase space, assuming the parent particle is at rest
 */

#ifndef GENERATOR
#define GENERATOR

#include<vector>
#include"TLorentzVector.h"
#include"TGenPhaseSpace.h"
 
class Generator {
  public:
    /**
     * Constructor that takes in the particle passes and sets up phase space
     * @param mass_parent Mass of parent particle
     * @param mass_decay mass of decay particles
     * @param particles Number of particles in the final state
     */
    Generator(const double &mass_parent, const Double_t *mass_decay, Int_t particles);
    /**
     * Function that generates a random unweighted event
     */
    std::vector<TLorentzVector> Generate();
 private:
    /**
     * TGenPhaseSpace object that describes phase space of decay particles
     */
    TGenPhaseSpace m_phasespace;
};

#endif
