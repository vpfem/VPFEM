#ifndef CONCRETE01
#define CONCRETE01

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cmath>

#include "ElasticPlainStress.h"
#include "Log/Log.h"

// This class models the reinforced concrete based on the paper 
//"reinforced concrete element memberane element formulation"
// Reinforced concrete model based on the MCFT method
// This class is designed for steel in x and y direction 
class Concrete01 : public ElasticPlainStress {
private:
    double fcPrime;         // fcPrime: f'c of concrete
    double fcr;             // fcr : tensile capacity of concrete
    double Ec;              // Ec : elastic modulus of concrete
    double e0;              // e0 : strain ar f'c
    double nu;              // nu : possion ratio
    double Es;              // Es : elastic modulus of steel
    double beta;             // beta: the angle of the element
    double ro_x;            // ro_x : steel ratio in x direction
    double fyx;             // fyx : yield stress of steel
    double alpha_x;         // aplha_x: angle of the rebars in x direction
    double ro_y;            // ro_y : steel ratio in y direction
    double fyy;             // fyy : yield stress of steel
    double alpha_y;         // aplha_y: angle of the rebars in y direction
    // cracking parameters
    double e_cr;
public:
    Concrete01(double, double, double, double, double, double, double, double, double, double, double, double, double);
    ~Concrete01();
    double UpdateMatMatrix(double*, double*, double*) override;
    void PrincipalStrain(double* result, double& theta, double* e);
    void SteelStrain(double* result, double* strain, double beta, double alpha_x, double alpha_y);
    double SteelStrain(double* strain, double beta, double alpha);
    void StrainToStressConcerete(double *result, double* strain_p);
    void StrainToStressSteel(double* result, double* strain_p);
    void matMatrixConcrete(double* mat, double* strain, double* stress);
    void matMatrixSteel(double* mat_x, double* mat_y, double* strain, double* stress);
    void transformation(double* mat, double angle);
    double errorCheck(double E1, double E2, double Es1, double Es2, double*);
};
#endif