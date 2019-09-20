#ifndef STIFFNESSMATRIXFIRSTORDER_H
#define STIFFNESSMATRIXFIRSTORDER_H

#include <iostream>
#include <cuda.h>

#include "../StiffnessMatrix.h"
#include "../../Log/Log.h"
#include "../../Sparse/Sparse.h"

class StiffnessMatrixFirstOrder : public StiffnessMatrix
{
    //VARIABLES

    //METHODS
public:
    StiffnessMatrixFirstOrder();
    StiffnessMatrixFirstOrder(double*, Geometry&,   int);
    virtual ~StiffnessMatrixFirstOrder();
    virtual Sparse& GetStiffnessMatrix() = 0;
    int GetStiffnessMatrixSize() override;
    CUDA_HOSTDEV void stiffnessMatrixCalculation(  int,   int, double*,   int*, double*, double*, double*,   int* ,double*,   int*,   int*,   int*, const double*);
    CUDA_HOSTDEV void constantCreator(  int numberElement, double* c, double* x, double* y,   int* mesh);
    CUDA_HOSTDEV void DisplacementToStrainEachElement(int, double *,   int *,   int *, double*, int,   int ,double* ,   int* , double* );
    virtual void DisplacementToStrain(double* displacement, double* strain) = 0;
};
#endif


