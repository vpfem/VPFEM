#pragma once 

#include <iostream>
#include "ConjugateGradient.h"
#include "ConjugateGradientPrecond.h"
#include "Cholesky.h"
#include "Log/Log.h"

class Solver {
    protected:
        static constexpr double tol = 0.00001;
        static constexpr double max_iter = 10000;
        int N, nz;
        double *val;
        int *rowPtr, *colIndex;
        double *x, *rhs, *rhs_copy;

    public:
        Solver(int N, int nz, double * val, int *rowPtr, int *colIndex, double * x, double* rhs);
        ~Solver();
        // Direct Cholesky Decomposition for Sparse CSR matrices on GPU
        void DirectCholDevice();
        void IterativeCgDevice();
        void IterativeCgHost();
        void IterativeCgPrecondCholDevice();
        void IterativeCgPrecondLuDevice();
    private:
};
