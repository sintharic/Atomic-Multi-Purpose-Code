#include "globals.h"

void potentialType::updateCutOff() {
  if(atomProperties.size()>1) terminate("Mixing rules not yet implemented");
  
  // scale cut off with current volumetric strain
  double relPrefac = relCutOff;
  relPrefac *= pow(det3x3(add3x3(eye3,strain)), 1./3);

  rCutOff2.clear();
  if (atomProperties[0].r0crystal > 0)
    rCutOff2.push_back( {pow(relPrefac * atomProperties[0].r0crystal, 2)} );
  else
    rCutOff2.push_back( {pow(relPrefac * atomProperties[0].r0dimer, 2)} );
};

void updateChargeDens(potentialType);
void mixingRuleC6(potentialType&);
void mixingRuleLJ(potentialType&);
void mixingRuleEAM(potentialType&);
void mixingRuleATM(potentialType&);

double energyC6(int, potentialType);
double forceC6(int, potentialType);
double energyLJ(int, potentialType);
double forceLJ(int, potentialType);
double energyEAM(int, potentialType);
double forceEAM(int, potentialType);
double energyATM(int, potentialType);
double forceATM(int, potentialType);