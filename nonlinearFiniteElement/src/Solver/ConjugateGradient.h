#pragma once 

#include <iostream>
#include <string>

/* Using updated (v2) interfaces to cublas */
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>

#include "Log/Log.h"

#include "SparseMatrixOperations.h"

namespace GPU {
    class ConjugateGradient {
    private:

    public:
        ConjugateGradient(int N, int nz, double * val, int *rowPtr, int *colIndex, double * x, double* rhs, const double tol);
        ~ConjugateGradient();
    private:
    };
}

namespace SingleCPU {
    class ConjugateGradient {
    private:

    public:
        ConjugateGradient(int N, int nz, double * val, int *rowPtr, int *colIndex, double * x, double* rhs, const double tol);
        ~ConjugateGradient();
    private:
    };
}