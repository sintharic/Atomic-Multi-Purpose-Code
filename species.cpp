#include "globals.h"

void initSpecies() {

  // Lennard-Jonesium
  species.LennardJonesium.name = "LJ";
  species.LennardJonesium.chemSymb = "Ar";
  species.LennardJonesium.mass = 1;
  species.LennardJonesium.charge = 0;
  species.LennardJonesium.scattAmpl = Complex(1,0);
  species.LennardJonesium.u0dimer = 1;
  species.LennardJonesium.r0dimer = 1;
  species.LennardJonesium.C6 = 1;
  species.LennardJonesium.C9 = 1;
  species.LennardJonesium.crystal = "fcc";
  species.LennardJonesium.LJ.u0 = species.LennardJonesium.u0dimer;
  species.LennardJonesium.LJ.r0 = species.LennardJonesium.r0dimer;

  // Ducastellium
  species.Ducastellium.name = "Ducastellium";
  species.Ducastellium.chemSymb = "Cu";
  species.Ducastellium.mass = 1;
  species.Ducastellium.charge = 0;
  species.Ducastellium.scattAmpl = Complex(1,0);
  species.Ducastellium.C6 = 1;
  species.Ducastellium.C9 = 1;
  species.Ducastellium.crystal = "fcc";
  species.Ducastellium.u0crystal = 1.0;
  species.Ducastellium.r0crystal = 1.0;
  species.Ducastellium.EAM.sigmaRep = 0.1;
  species.Ducastellium.EAM.sigmaAtt = 0.25;
  species.Ducastellium.EAM.vRep = 918;
  species.Ducastellium.EAM.aAtt = 7.11;

  // argon
  species.argon.name = "Argon";
  species.argon.chemSymb = "Ar";
  species.argon.mass = 39.948;
  species.argon.charge = 0;
  species.argon.polarizability = 0;
  species.argon.scattAmpl = Complex(1,0);
  species.argon.u0dimer = 120*kBoltzK;
  species.argon.r0dimer = 3.5;
  species.argon.LJ.u0 = species.argon.u0dimer;
  species.argon.LJ.r0 = species.argon.r0dimer;
  getC6fromLJ(species.argon);
  //species.argon.C6 = 2*species.argon.u0dimer*pow(species.argon.r0dimer,6);
  species.argon.crystal = "fcc";
  species.argon.r0crystal = 3.3;
};

void getC6fromLJ(speciesType &spec) {
  spec.C6 = 2 * spec.LJ.u0 * pow(spec.LJ.r0, 6);
}

void getC9fromC6(speciesType &spec) {
  if (spec.C6 == 1)
    terminate("Cannot reassign C9 coefficient for generic atom type.");
  spec.C9 = spec.C6 * sqrt(spec.C6/hartree);
}