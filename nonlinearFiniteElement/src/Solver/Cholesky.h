#pragma once 

#include <iostream>

/* Using updated (v2) interfaces to cublas */
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cusolverSp.h>

#include "Log/Log.h"
namespace GPU 
{
    class Cholesky  
    {
    private:


    public:
      Cholesky(int N, int nz, double * val, int *rowPtr, int *colIndex, double * x, double* rhs, const double tol);
      ~Cholesky();
    private:

    };
}