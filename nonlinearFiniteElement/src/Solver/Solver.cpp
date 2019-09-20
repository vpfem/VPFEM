#include "Solver.h"

Solver::Solver(int N, int nz, double * val, int *rowPtr, 
                int *colIndex, double * x, double* rhs) 
: N(N), nz(nz), val(val), rowPtr(rowPtr), 
colIndex(colIndex), x(x), rhs(rhs)
{
    cudaMallocManaged(&rhs_copy,N*sizeof(double));
    cudaMemcpy(rhs_copy, rhs, N*sizeof(double),cudaMemcpyDeviceToDevice);
}

/*

    INFO("Conjugate Gradient Precondition Selected (LU) (GPU)");
    cudaMemset(x, 0, N*sizeof(double));
    cudaMemcpy(rhs, rhs, N*sizeof(double),cudaMemcpyDeviceToDevice);
    GPU::ConjugateGradientPrecondLU(N, nz, val, rowPtr, colIndex, x, rhs, tol);
    cudaFree(rhs);
    INFO("Conjugate Gradient Precondition Selected (CHOL) (GPU)");
    cudaMemset(x, 0, N*sizeof(double));
    //cudaMemcpy(rhs, rhs, N*sizeof(double),cudaMemcpyDeviceToDevice);
    GPU::ConjugateGradientPrecondChol(N, nz, val, rowPtr, colIndex, x, rhs, tol);
    cudaFree(rhs);
*/
void Solver::DirectCholDevice() 
{
    INFO("Cholesky Decomposition Selected");
    cudaMemset(x, 0, N*sizeof(double));
    GPU::Cholesky(N, nz, val, rowPtr, colIndex, x, rhs_copy, tol);
}

void Solver::IterativeCgDevice()
{
    INFO("Conjugate Gradient Selected (GPU)");
    cudaMemset(x, 0, N*sizeof(double));
    //ARRAY(rhs,N);
    GPU::ConjugateGradient(N, nz, val, rowPtr, colIndex, x, rhs_copy, tol);
}

void Solver::IterativeCgHost()
{
    INFO("Conjugate Gradient Selected (SingleCPU)");
    cudaMemset(x, 0, N*sizeof(double));
    SingleCPU::ConjugateGradient(N, nz, val, rowPtr, colIndex, x, rhs_copy, tol);
}

void Solver::IterativeCgPrecondCholDevice()
{
    INFO("Conjugate Gradient with LU Preconditioning Selected (GPU)");
    cudaMemset(x, 0, N*sizeof(double));
    GPU::ConjugateGradientPrecondChol(N, nz, val, rowPtr, colIndex, x, rhs_copy, tol);
}
void Solver::IterativeCgPrecondLuDevice()
{
    INFO("Conjugate Gradient with LU Preconditioning Selected (GPU)");
    cudaMemset(x, 0, N*sizeof(double));
    GPU::ConjugateGradientPrecondLU(N, nz, val, rowPtr, colIndex, x, rhs_copy, tol);
}

Solver::~Solver() {
}