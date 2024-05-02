#ifndef HEADER_H
#define HEADER_H
//endif at EOF

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <complex>
#include <string>
#include <limits.h>

using namespace std;

#include "unit.h"

const double INF = std::numeric_limits<double>::infinity();


typedef ptrdiff_t Lint;
typedef std::complex<double> Complex;

template <class T> int sign(T a){
  return(a<0 ? -1 : 1);
}

void terminate(string msg) {
  cerr << msg << endl;
  exit(1);
}

typedef array<double, 3> Vec3D;
typedef vector<Vec3D> MatrixNx3;



// ----- Vec3D type for position vectors -----

Vec3D operator+(Vec3D v1, Vec3D v2) {
  Vec3D result;
  for (int iDim = 0; iDim < 3; ++iDim) result[iDim] = v1[iDim] + v2[iDim];
  return result;
}

Vec3D operator-(Vec3D v1, Vec3D v2) {
  Vec3D result;
  for (int iDim = 0; iDim < 3; ++iDim) result[iDim] = v1[iDim] - v2[iDim];
  return result;
}

Vec3D operator*(double val, Vec3D vec) {
  Vec3D result;
  for (int iDim = 0; iDim < 3; ++iDim) result[iDim] = val*vec[iDim];
  return result;
}

Vec3D operator/(Vec3D vec, double val) {
  Vec3D result;
  for (int iDim = 0; iDim < 3; ++iDim) result[iDim] = vec[iDim]/val;
  return result;
}

double operator*(Vec3D v1, Vec3D v2) { // scalar product
  return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}

MatrixNx3 operator+(MatrixNx3 Mat, Vec3D vec) {
  MatrixNx3 result = Mat;
  for (int iRow = 0; iRow < Mat.size(); ++iRow) {
    Mat[iRow] = Mat[iRow] + vec;
  }
  return result;
}
MatrixNx3 operator+(Vec3D vec, MatrixNx3 Mat) {
  return Mat + vec;
}

MatrixNx3 operator*(double v, MatrixNx3 Mat) {
  MatrixNx3 result = Mat;
  for (int iRow = 0; iRow < Mat.size(); ++iRow) {
    for (int jCol = 0; jCol < 3; ++jCol) {
      result[iRow][jCol] *= v;
    }
  }
  return result;
}

// multiply a 3x3 matrix with a 3-element vector
Vec3D mult3x3Vec3D(MatrixNx3 Mat, Vec3D vec) {
  if (Mat.size() != 3) terminate("Matrix with shape != 3x3 encountered in mult3x3Vec3D().");

  Vec3D result;
  for (int iDim = 0; iDim < 3; ++iDim) {
    result[iDim] = Mat[iDim][0]*vec[0] + Mat[iDim][1]*vec[1] + Mat[iDim][2]*vec[2];
  }

  return result;
}

MatrixNx3 mult3x3(MatrixNx3 M1, MatrixNx3 M2) {
  if (M1.size() != 3 || M2.size() != 3) terminate("Matrix with shape != 3x3 encountered in mult3x3().");

  MatrixNx3 M3(3);
  for (int iDim = 0; iDim < 3; ++iDim) {
    for (int jDim = 0; jDim < 3; ++jDim) {
      for (int kDim = 0; kDim < 3; ++kDim) {
        M3[iDim][jDim] += M1[iDim][kDim]*M2[kDim][jDim];
      }
    }
  }

  return M3;
}

MatrixNx3 add3x3(MatrixNx3 M1, MatrixNx3 M2) { // scalar product
  if (M1.size() != M2.size()) terminate("Cannot add matrices shapes "+to_string(M1.size())+"x3 and "+to_string(M2.size())+"x3.");

  MatrixNx3 M3(M1.size());
  for (int iRow = 0; iRow < M1.size(); ++iRow) {
    for (int jCol = 0; jCol < 3; ++jCol) {
      M3[iRow][jCol] = M1[iRow][jCol] + M2[iRow][jCol];
    }
  }

  return M3;
}

double det3x3(MatrixNx3 Mat) { // scalar product
  if (Mat.size() != 3) terminate("Matrix with shape != 3x3 encountered in det3x3().");

  double result = Mat[0][0] * (Mat[1][1] * Mat[2][2] - Mat[2][1] * Mat[1][2]) -
                  Mat[0][1] * (Mat[1][0] * Mat[2][2] - Mat[1][2] * Mat[2][0]) +
                  Mat[0][2] * (Mat[1][0] * Mat[2][1] - Mat[1][1] * Mat[2][0]);
  return result;
}


// ----- structure for interaction parameters -----

struct paramsLJ {
  double u0, r0;
};

struct paramsBuckingham {
  double vRep, sigmaRep;
};

struct paramsEAM {
  double vRep, sigmaRep, aAtt, sigmaAtt, cd;
};

struct paramsQE {
  double hard, EN;
};



// ----- structure for species type (atom, ion) parameters -----

struct speciesType {
  string name, chemSymb;

  // physical properties
  double mass = 1, charge = 0, polarizability = 0;
  Complex scattAmpl = Complex(0,1);
  double u0dimer = 0, r0dimer = 0;
  double C6 = 0, C9 = 0; // can be chosen different from dimer values

  // equilibrium crystal properties
  string crystal = "none";
  double u0crystal = 0, r0crystal = 0;

  // interatomic interactions
  paramsLJ LJ;
  paramsBuckingham Buckingham;
  paramsEAM EAM;
  paramsQE QE;
};



// ----- structure for crystal type -----

struct crystalType {
  string name;
  Vec3D latticeConstant;
  MatrixNx3 basis;
  vector<int> atomType;
};



// ----- class for potential type -----

class potentialType {
public:
  string name;

  // whether or not this potential requires chargeDensities
  bool fUpdateCD = 0;

  // potential parameters
  vector< vector<double> > C6;
  vector< vector<double> > u0, r0sq;
  vector< vector<double> > vRep, rangeRep;
  vector< vector<double> > aAtt, vAtt, rangeAtt;
  vector< vector<double> > C9;

  // cut-off parameters
  double relCutOff = INF;
  vector< vector<double> > rCutOff2;
  void updateCutOff();

  // main functions
  double (*energyFunction)(int, potentialType);
  double energy(int iAtom) {return energyFunction(iAtom, (*this));}
  double (*forceFunction)(int, potentialType);
  double force(int iAtom) {return forceFunction(iAtom, (*this));}
  void (*mixingFunction)(potentialType&);
  void mixingRule() {mixingFunction((*this));}

  //potentialType(double (*energy)(int), double (*force)(int)) {
  //  (*this).energy = energy;
  //  (*this).force = force;
  //};
  //static double energyC6(int);
  //static double forceC6(int);
  //static double energyLJ(int);
  //static double forceLJ(int);
  //static double energyEAM(int);
  //static double forceEAM(int);
  //static double energyATM(int);
  //static double forceATM(int);
};

#endif