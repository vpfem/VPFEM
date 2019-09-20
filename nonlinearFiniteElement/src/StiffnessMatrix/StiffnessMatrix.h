#ifndef STIFFNESSMATRIX_H
#define STIFFNESSMATRIX_H

//macro
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

#include <iostream>
#include <string>
#include <cuda_runtime_api.h>
#include <cuda.h>

#include "../Material/Material.h"
#include "../Geometry/Geometry.h"
#include "../Log/Log.h"
#include "../Timer/Timer.h"
#include "../Sparse/Sparse.h"
#include "SparseMatrix/SparseMatrix.h"

class StiffnessMatrix
{
    //VARIABLES
public:
    Sparse* stiffMat;
protected:
    static const   int dimention = 2;
      int numberOfIntegrationPoint;
      int numberOfElements; //number of elements
      int nipSquared; // number of interation point squared
    double* integrationNode;   int* integrationPos; double* integrationWeight;
    double* material;
    Geometry* geometry;
    double* c;
      int sizeStiffMatPerEle;
      int simulationSize;
      int stiffMatSize;
    //std::stirng hardwareType;
private:
    //METHODS
public:
    StiffnessMatrix(double*, Geometry&,   int);
    virtual ~StiffnessMatrix();
    virtual int GetStiffnessMatrixSize() = 0;
    virtual Sparse& GetStiffnessMatrix() = 0;

protected:
private:
    CUDA_HOSTDEV void integrationPoint();
};
#endif
