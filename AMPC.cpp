// atomic multi-purpose code

#include "AMPC.h"

int main(){
  initParams();

  cout << indivEnergy(0) << endl;
  //findLatticeConst();
  //elasticTensor();
}

void initParams(){
  initSpecies();
  initCrystals();
  initPotentials();

  // define atom species: choose pre-defined (and alter it) or make a new one 
  speciesType s = species.LennardJonesium;
  atomProperties.push_back(s);

  // create crystal: choose pre-defined (and alter it) or make a new one 
  crystalType c = crystals.fcc;
  makeCrystal(c,36,36,36);

  // define interactions: choose pre-defined (and alter it) or make a new one 
  potentialType p = potentials.C6;
  p.relCutOff = INF;
  interactions.push_back(p);

}

void findLatticeConst() {
  MatrixNx3 hMat0 = hMat;
  double delta_eps = 0.02;
  cout << "tr(epsilon)\tU_c" << endl;
  for (double tr_eps = -20*delta_eps; tr_eps <= 20.001*delta_eps; tr_eps+=delta_eps) {
    for (int iDim = 0; iDim < nDim; ++iDim) {
      strain[0][0] = tr_eps/3;
      strain[1][1] = tr_eps/3;
      strain[2][2] = tr_eps/3;
      hMat[iDim][iDim] = hMat0[iDim][iDim]*(1+tr_eps/3);
      makeMetric();
    }
    cout << tr_eps << "\t" << totalEnergy()/nAtom << endl;
  }
}

void elasticTensor() {
  double delta_epsilon = 0.002;
  MatrixNx3 hMat0 = hMat;
  double Volume0 = det3x3(hMat);
  for (int iDim = 0; iDim < nDim; ++iDim) {
    for (int jDim = iDim; jDim < nDim; ++jDim) {
      vector<double> FreeEnergy = {};
      for (double eps_ij = -delta_epsilon; eps_ij <= delta_epsilon; eps_ij+=delta_epsilon) {
        strain[iDim][jDim] = eps_ij;
        strain[jDim][iDim] = eps_ij;
        hMat = mult3x3(add3x3(eye3,strain),hMat0);
        makeMetric();
        double Volume = det3x3(hMat);
        FreeEnergy.push_back(totalEnergy()/Volume0);
      }
      printf("C_%i%i%i%i = %.8g\n", iDim+1, jDim+1, iDim+1, jDim+1, 
        (FreeEnergy[2] + FreeEnergy[0] - 2*FreeEnergy[1]) / (delta_epsilon*delta_epsilon));
      strain[iDim][jDim] = 0;
      strain[jDim][iDim] = 0;
    }
  }
}
