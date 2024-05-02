#ifndef GLOBALS_H
#define GLOBALS_H
//endif at EOF

#include "header.h"

// mini database of atom/ion species (see species.cpp)
struct speciesList {
  speciesType LennardJonesium, Ducastellium;
  speciesType argon, copper, gold, aluminum;
};
speciesList species;

// mini database of crystal structures (see crystals.cpp)
struct crystalList {
  crystalType sc, bcc, fcc, hcp, dia;
  crystalType NaCl, CsCl;
};
crystalList crystals;

// mini database of potentials (see interactions.cpp)
struct potentialList {
  potentialType C6;
  potentialType LJ;
  potentialType EAM;
  potentialType ATM;
};
potentialList potentials;



// fields of size nAtom
int nAtom;
vector<int> atomType;
vector<double> atomCharge, chargeDensity;
MatrixNx3 atomFracPos, atomRealPos;

// properties of all included atomTypes
vector <speciesType> atomProperties;



// system dimensions (config.cpp)
const int nDim = 3;
MatrixNx3 hMat(nDim), metric(nDim);
MatrixNx3 eye3 = {{1,0,0},{0,1,0},{0,0,1}};
MatrixNx3 strain = MatrixNx3(3);
inline Vec3D PBC(Vec3D);
double getSqDist(int,int,bool);
void trans2Real(int);
void makeMetric();



// functions for specieTypes (species.cpp)
void initSpecies();
void getC6fromLJ(speciesType&);
void getC9fromC6(speciesType&);

// functions for crystalTypes (crystal.cpp)
void initCrystals();
void makeCrystal(crystalType crys, int, int, int);
void dumpCrystal();

// functions for potentialTypes (interaction.cpp)
vector <potentialType> interactions;
void initPotentials();
double totalEnergy();
double indivEnergy(int);
double allForces();




// apply periodic boundary conditions to a fractional position vector
inline Vec3D PBC(Vec3D vec) {
  Vec3D result = vec;
  for (int iDim = 0; iDim < nDim; ++iDim) {
    if (result[iDim] < -0.5) result[iDim] += 1;
    else if (result[iDim] > 0.5) result[iDim] -= 1;
  }
  return result;
}

// getSqDist is efficient if the actual vector r_ij is not of interest
double getSqDist(int iAtom, int jAtom, bool fShear){
  Vec3D dr = PBC(atomFracPos[iAtom] - atomFracPos[jAtom]);

  double r2 = 0;
  for (int iDim=0; iDim<nDim; ++iDim)
    r2 += metric[iDim][iDim]*dr[iDim]*dr[iDim];

  if (!fShear) return(r2);

  r2 += metric[0][1]*dr[0]*dr[1];
  r2 += metric[0][2]*dr[0]*dr[2];
  r2 += metric[1][2]*dr[1]*dr[2];

  return(r2);
}

void trans2Real(int iAtom){
  atomRealPos[iAtom] = mult3x3Vec3D(hMat, atomFracPos[iAtom]);
}

#endif