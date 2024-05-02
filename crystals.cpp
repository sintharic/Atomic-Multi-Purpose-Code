#include "globals.h"


// compute metric tensor for more efficient squares of vectors
void makeMetric(){  // constructs metric tensor, upper diagonal only
  for (int iDim=0; iDim<nDim; ++iDim) {
  for (int jDim=iDim; jDim<nDim; ++jDim) {
    metric[iDim][jDim] = metric[jDim][iDim] = 0;
    for (int kDim=0; kDim<nDim; ++kDim) {
      metric[iDim][jDim] += hMat[kDim][iDim]*hMat[kDim][jDim];
    }
    if (iDim!=jDim) metric[iDim][jDim] *= 2;
  } }
}



void makeCrystal(crystalType crys, int nx, int ny, int nz){

  double unitLength = 1;
  if (atomProperties[0].r0crystal > 0) unitLength = atomProperties[0].r0crystal;
  else if (atomProperties[0].r0dimer > 0) unitLength = atomProperties[0].r0dimer;
  else terminate("No length unit defined in species " + atomProperties[0].name + ".");


  // fill h matrix and metric tensor
  hMat[0][0] = unitLength*crys.latticeConstant[0]*nx;
  hMat[1][1] = unitLength*crys.latticeConstant[1]*ny;
  hMat[2][2] = unitLength*crys.latticeConstant[2]*nz;
  makeMetric();

  // reset atom positions and properties
  nAtom = nx*ny*nz*crys.basis.size();
  atomFracPos.clear();
  atomRealPos.clear();
  atomCharge.clear();
  atomType.clear();

  // fill atom positions and properties
  for (int ix = 0; ix < nx; ++ix) {
  for (int iy = 0; iy < ny; ++iy) {
  for (int iz = 0; iz < nz; ++iz) {
    Vec3D posCell = {1.*ix, 1.*iy, 1.*iz};
    for (int iB = 0; iB < crys.basis.size(); ++iB) {
      Vec3D pos = posCell + crys.basis[iB];
      pos[0] /= nx;
      pos[1] /= ny;
      pos[2] /= nz;
      atomFracPos.push_back(pos);
      atomCharge.push_back(atomProperties[crys.atomType[iB]].charge);
      atomType.push_back(crys.atomType[iB]);
      chargeDensity.push_back(0.);
    }
  } } }

  atomRealPos.resize(nAtom);
  for (int iAtom=0; iAtom<nAtom; ++iAtom) trans2Real(iAtom);

}


void initCrystals() {
  Vec3D e1 = {1, 0, 0}, e2 = {0, 1, 0}, e3 = {0, 0, 1};

  // simple cubic
  crystals.sc.name = "sc";
  crystals.sc.latticeConstant[0] = 1.;
  crystals.sc.latticeConstant[1] = 1.;
  crystals.sc.latticeConstant[2] = 1.;
  crystals.sc.basis.push_back(0.*e1);
  crystals.sc.atomType.assign(crystals.sc.basis.size(), 0);

  // body-centered cubic
  crystals.bcc.name = "bcc";
  crystals.bcc.latticeConstant[0] = 2./sqrt(3.);
  crystals.bcc.latticeConstant[1] = 2./sqrt(3.);
  crystals.bcc.latticeConstant[2] = 2./sqrt(3.);
  crystals.bcc.basis.push_back(0. *e1 + 0. *e2 + 0. *e3);
  crystals.bcc.basis.push_back(0.5*e1 + 0.5*e2 + 0.5*e3); 
  crystals.bcc.atomType.assign(crystals.bcc.basis.size(), 0);

  // face-centered cubic
  crystals.fcc.name = "fcc";
  crystals.fcc.latticeConstant[0] = sqrt(2.);
  crystals.fcc.latticeConstant[1] = sqrt(2.);
  crystals.fcc.latticeConstant[2] = sqrt(2.);
  crystals.fcc.basis.push_back(0. *e1 + 0. *e2 + 0. *e3);
  crystals.fcc.basis.push_back(0.5*e1 + 0.5*e2 + 0. *e3);
  crystals.fcc.basis.push_back(0.5*e1 + 0. *e2 + 0.5*e3);
  crystals.fcc.basis.push_back(0. *e1 + 0.5*e2 + 0.5*e3);
  crystals.fcc.atomType.assign(crystals.fcc.basis.size(), 0);

  // hexagonal close-packed
  crystals.hcp.name = "hcp";
  crystals.hcp.latticeConstant[0] = 1.;
  crystals.hcp.latticeConstant[1] = sqrt(3.);
  crystals.hcp.latticeConstant[2] = 2*sqrt(2./3);
  crystals.hcp.basis.push_back(0. *e1 + 0.  *e2 + 0. *e3);
  crystals.hcp.basis.push_back(0.5*e1 + 0.25*e2 + 0.5*e3);
  crystals.hcp.basis.push_back(1.0*e1 + 0.75*e2 + 0.5*e3);
  crystals.hcp.basis.push_back(0.5*e1 + 0.5 *e2 + 0. *e3);
  crystals.hcp.atomType.assign(crystals.hcp.basis.size(), 0);

  // diamond
  crystals.dia.name = "diamond";
  crystals.dia.latticeConstant[0] = 4./sqrt(3.);
  crystals.dia.latticeConstant[1] = 4./sqrt(3.);
  crystals.dia.latticeConstant[2] = 4./sqrt(3.);
  crystals.dia.basis.push_back(0.  *e1 + 0.  *e2 + 0.  *e3);
  crystals.dia.basis.push_back(0.5 *e1 + 0.5 *e2 + 0.  *e3);
  crystals.dia.basis.push_back(0.5 *e1 + 0.  *e2 + 0.5 *e3);
  crystals.dia.basis.push_back(0.  *e1 + 0.5 *e2 + 0.5 *e3);
  crystals.dia.basis.push_back(0.25*e1 + 0.25*e2 + 0.25*e3);
  crystals.dia.basis.push_back(0.25*e1 + 0.75*e2 + 0.75*e3);
  crystals.dia.basis.push_back(0.75*e1 + 0.25*e2 + 0.75*e3);
  crystals.dia.basis.push_back(0.75*e1 + 0.75*e2 + 0.25*e3);
  crystals.dia.atomType.assign(crystals.dia.basis.size(), 0);

  // CsCl
  crystals.CsCl = crystals.bcc;
  crystals.CsCl.name = "CsCl";
  crystals.CsCl.latticeConstant[0] = 2./sqrt(3.);
  crystals.CsCl.latticeConstant[1] = 2./sqrt(3.);
  crystals.CsCl.latticeConstant[2] = 2./sqrt(3.);
  crystals.CsCl.atomType[1] = 1;

  // NaCl
  crystals.NaCl = crystals.fcc;
  crystals.NaCl.name = "NaCl";
  crystals.NaCl.latticeConstant[0] = 2.;
  crystals.NaCl.latticeConstant[1] = 2.;
  crystals.NaCl.latticeConstant[2] = 2.;
  crystals.NaCl.basis.push_back(0.5*e1 + 0. *e2 + 0. *e3);
  crystals.NaCl.basis.push_back(0. *e1 + 0.5*e2 + 0. *e3);
  crystals.NaCl.basis.push_back(0. *e1 + 0. *e2 + 0.5*e3);
  crystals.NaCl.basis.push_back(0.5*e1 + 0.5*e2 + 0.5*e3);
  for (int i = 0; i < 4; ++i) crystals.NaCl.atomType.push_back(1);

}





// old initialization routines for reference


// ----- structure for crystal type -----
/*
struct crystalType {
  string name;
  Vec3D latticeConstant;
  MatrixNx3 basis;
  vector<int> atomType;
};

crystalType CsCl = {
  "CsCl", 
  {2./sqrt(3.), 2./sqrt(3.), 2./sqrt(3.)},
  {{0  , 0  , 0  },
   {0.5, 0.5, 0.5}},
  {0,1}
};

crystalType NaCl = {
  "NaCl", 
  {2., 2., 2.},
  {{0  , 0  , 0  },
   {0.5, 0.5, 0  },
   {0.5, 0  , 0.5},
   {0  , 0.5, 0.5},
   {0.5, 0  , 0  },
   {0. , 0.5, 0  },
   {0  , 0  , 0.5},
   {0.5, 0.5, 0.5}},
   {0,0,0,0,1,1,1,1}
};

crystalType sc = {
  "sc", 
  {1., 1., 1.},
  {{0, 0, 0}},
  {0},
};

crystalType bcc = {
  "bcc", 
  {2./sqrt(3.), 2./sqrt(3.), 2./sqrt(3.)},
  {{0  , 0  , 0  },
   {0.5, 0.5, 0.5}},
  {0,0},
};

crystalType fcc = {
  "fcc", 
  {sqrt(2.), sqrt(2.), sqrt(2.)},
  {{0  , 0  , 0  },
   {0.5, 0.5, 0  },
   {0.5, 0  , 0.5},
   {0  , 0.5, 0.5}},
  {0,0,0,0},
};

crystalType hcp = {
  "hcp", 
  {1, sqrt(3.), 2*sqrt(2./3)},
  {{0    , 0    , 0    },
   {0.5  , 0.25 , 0.5  },
   {1    , 0.75 , 0.5  },
   {0.5  , 0.5  , 0    }},
  {0,0,0,0},
};

crystalType dia = {
  "diamond", 
  {4./sqrt(3.), 4./sqrt(3.), 4./sqrt(3.)},
  {{0   , 0   , 0   },
   {0.5 , 0.5 , 0   },
   {0.5 , 0   , 0.5 },
   {0   , 0.5 , 0.5 },
   {0.25, 0.25, 0.25},
   {0.25, 0.75, 0.75},
   {0.75, 0.25, 0.75},
   {0.75, 0.75, 0.25}},
  {0,0,0,0,0,0,0,0},
};
*/