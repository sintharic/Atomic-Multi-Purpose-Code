#include "interactions.h"


void initPotentials() {

  potentials.C6.name = "C6";
  potentials.C6.energyFunction = energyC6;
  potentials.C6.forceFunction = forceC6;
  potentials.C6.mixingFunction = mixingRuleC6;

  potentials.LJ.name = "LJ";
  potentials.LJ.energyFunction = energyLJ;
  potentials.LJ.forceFunction = forceLJ;
  potentials.LJ.mixingFunction = mixingRuleLJ;

  potentials.EAM.name = "EAM";
  potentials.EAM.energyFunction = energyEAM;
  potentials.EAM.forceFunction = forceEAM;
  potentials.EAM.mixingFunction = mixingRuleEAM;
  potentials.EAM.relCutOff = 1.1;
  potentials.EAM.fUpdateCD = 1;

  potentials.ATM.name = "ATM";
  potentials.ATM.energyFunction = energyATM;
  potentials.ATM.forceFunction = forceATM;
  potentials.ATM.mixingFunction = mixingRuleATM;

}


double indivEnergy(const int iAtom) {
  double potEnergy = 0;
  for (potentialType inter : interactions) {
    inter.updateCutOff();
    inter.mixingRule();
    if (inter.fUpdateCD) updateChargeDens(inter);
    potEnergy += inter.energy(iAtom);
  }
  return potEnergy;
}

double totalEnergy(){
  double potEnergy = 0;
  for (potentialType inter : interactions) {
    inter.updateCutOff();
    inter.mixingRule();
    if (inter.fUpdateCD) updateChargeDens(inter);
    for (int iAtom=0; iAtom<as_const(nAtom); ++iAtom) {
      potEnergy += inter.energy(iAtom);
    }
  }
  return(potEnergy);
}

double allForces(){
  double potEnergy = 0;
  for (potentialType inter : interactions) {
    for (int iAtom=0; iAtom<as_const(nAtom); ++iAtom) 
    potEnergy += inter.force(iAtom);
  }
  return(potEnergy);
}

double energyC6(const int iAtom, potentialType pot){

  const int iType = atomType[iAtom];

  double sum = 0;
  for (int jAtom=0; jAtom<as_const(nAtom); ++jAtom) {
    if (jAtom==iAtom) continue;

    int jType = atomType[jAtom];

    double rij2 = getSqDist(iAtom,jAtom,false);
    if (rij2>pot.rCutOff2[iType][jType]) continue;

    sum += pot.C6[iType][jType]/(rij2*rij2*rij2);
  }
  return(sum);
}

double forceC6(const int iAtom, potentialType pot){
  return(0);
}

void mixingRuleC6(potentialType &pot) {
  if (atomProperties.size() > 1) terminate("Mixing rules not yet implemented");
  pot.C6.clear();
  pot.C6.push_back({atomProperties[0].C6});
}

double energyLJ(const int iAtom, potentialType pot){

  const int iType = atomType[iAtom];

  double potEnergy = 0;
  for (int jAtom=iAtom+1; jAtom<as_const(nAtom); ++jAtom) {

    int jType = atomType[jAtom];
    double r2 = getSqDist(iAtom,jAtom,true);
    if (r2>pot.rCutOff2[iType][jType]) continue;

    r2 /= pot.r0sq[iType][jType];
    double r6Inv = 1./(r2*r2*r2);
    potEnergy += pot.u0[iType][jType]*(r6Inv*r6Inv-2*r6Inv);

  } 

  return(potEnergy);
}

double forceLJ(const int iAtom, potentialType pot){
  return(0);
}

void mixingRuleLJ(potentialType &pot) {
  if (atomProperties.size() > 1) terminate("Mixing rules not yet implemented");
  pot.u0.clear();
  pot.u0.push_back({atomProperties[0].LJ.u0});
  pot.r0sq.clear();
  pot.r0sq.push_back({pow(atomProperties[0].LJ.r0,2)});
}

void updateChargeDens(potentialType pot) {

  for (int iAtom=0; iAtom < as_const(nAtom); ++iAtom) chargeDensity[iAtom] = 0;

  for (int iAtom = 0; iAtom < as_const(nAtom); ++iAtom) {
    int iType = atomType[iAtom];
    for (int jAtom = 0; jAtom < nAtom; ++jAtom) {
      if (iAtom == jAtom) continue;
      int jType = atomType[jAtom];
      double rij2 = getSqDist(iAtom,jAtom,true);
      if (rij2 > pot.rCutOff2[iType][jType]) continue;

      chargeDensity[iAtom] += pot.aAtt[iType][jType]*exp(-sqrt(rij2)/pot.rangeAtt[iType][jType]);
    }
  }
}

double energyEAM(const int iAtom, potentialType pot) {

  double potEnergy = - sqrt(chargeDensity[iAtom]);
  int iType = atomType[iAtom];
  for (int jAtom = 0; jAtom < as_const(nAtom); ++jAtom) {
    if (iAtom==jAtom) continue;
    int jType = atomType[jAtom];
    double rij2 = getSqDist(iAtom,jAtom,true);
    if (rij2 > pot.rCutOff2[iType][jType]) continue;
    potEnergy += pot.vRep[iType][jType]*exp(-sqrt(rij2)/pot.rangeRep[iType][jType])/2;
  }
  return potEnergy;

}

double forceEAM(const int iAtom, potentialType pot) {
  return 0;
}

void mixingRuleEAM(potentialType &pot) {
  if (atomProperties.size() > 1) terminate("Mixing rules not yet implemented");
  pot.vRep.clear();
  pot.vRep.push_back({atomProperties[0].EAM.vRep});
  pot.rangeRep.clear();
  pot.rangeRep.push_back({atomProperties[0].EAM.sigmaRep});
  pot.aAtt.clear();
  pot.aAtt.push_back({atomProperties[0].EAM.aAtt});
  pot.rangeAtt.clear();
  pot.rangeAtt.push_back({atomProperties[0].EAM.sigmaAtt});
}

double energyATM(const int iAtom, potentialType pot) {

  const Vec3D iPos = atomFracPos[iAtom];
  const int iType = atomType[iAtom];

  double atmSum = 0;
  for (int jAtom=0; jAtom<nAtom-1; ++jAtom) {
    if (jAtom==iAtom) continue;

    int jType = atomType[jAtom];

    Vec3D jPos = atomFracPos[jAtom];
    Vec3D Rij = mult3x3Vec3D(hMat, PBC(jPos-iPos));
    double rij2 = Rij*Rij;
    if (rij2>pot.rCutOff2[iType][jType]) continue;

    for (int kAtom=jAtom+1; kAtom<nAtom  ; ++kAtom) {
      if (kAtom==iAtom) continue;

      int kType = atomType[kAtom];

      Vec3D kPos = atomFracPos[kAtom];
      Vec3D Rjk = mult3x3Vec3D(hMat, PBC(kPos-jPos));
      double rjk2 = Rjk*Rjk;
      if (rjk2>pot.rCutOff2[jType][kType]) continue;

      Vec3D Rik = Rij+Rjk;
      double rik2 = Rik*Rik;
      if (rik2>pot.rCutOff2[iType][kType]) continue;

      double rij2rjk2rik2 = rij2*rjk2*rik2;
      double scalI = Rij*Rik;
      double scalJ = Rjk*Rij;
      double scalK = Rik*Rjk;
      double scalTotEff = 1-3*scalI*scalJ*scalK/rij2rjk2rik2;

      atmSum += pot.C9[iType][jType]*scalTotEff/(rij2rjk2rik2*sqrt(rij2rjk2rik2));
    }
  }

  return(atmSum);
}

double forceATM(const int iAtom, potentialType pot) {
  return 0;
}

void mixingRuleATM(potentialType &pot) {
  if (atomProperties.size() > 1) terminate("Mixing rules not yet implemented");
  pot.C9.clear();
  pot.C9.push_back({atomProperties[0].C9});
}