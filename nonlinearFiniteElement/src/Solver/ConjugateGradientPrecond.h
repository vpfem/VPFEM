#pragma once 

#include <iostream>
#include <string>

/* Using updated (v2) interfaces to cublas */
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <helper_functions.h>  // helper for shared functions common to CUDA Samples
#include <helper_cuda.h>       // helper function CUDA error checking and initialization

#include "Log/Log.h"

namespace GPU {
    class ConjugateGradientPrecondLU {
    private:

    public:
        ConjugateGradientPrecondLU(int N, int nz, double * val, int *rowPtr, int *colIndex, double * x, double* rhs, const double tol);
        ~ConjugateGradientPrecondLU();

    };

    class ConjugateGradientPrecondChol {
    private:

    public:
        ConjugateGradientPrecondChol(int N, int nz, double * val, int *rowPtr, int *colIndex, double * x, double* rhs, const double tol);
        ~ConjugateGradientPrecondChol();

    };
}