#include "Cholesky.h"
namespace GPU {
    Cholesky::Cholesky (int N, int nz, double * val, int *rowPtr, int *colIndex, double * x, double* rhs, const double tol)
    {
        // --- Start the cuda sparse
        cusparseHandle_t handle; cusparseCreate(&handle);
        // --- Descriptor for sparse matrix A
        cusparseMatDescr_t descrA;      
        cusparseCreateMatDescr(&descrA);
        cusparseSetMatType      (descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase (descrA, CUSPARSE_INDEX_BASE_ZERO);
        // --- CUDA solver initialization
        cusolverSpHandle_t solver_handle;
        cusolverSpCreate(&solver_handle);
        //
        int reorder = 0; 
        int singularity = 0;
        cudaDeviceSynchronize();
        cusolverSpDcsrlsvcholHost(solver_handle, N, nz, descrA, val, rowPtr, colIndex, rhs, tol, reorder, x, &singularity);
        if (singularity != -1) {
            WARNING("[Warning]: Matrix was not positive definitive trying other solvers...");
        }
        //for (int i = 0; i < N; i++) {std::cout << "F[" << i << "]= " << x[i] << std::endl;}
    }

    Cholesky::~Cholesky()
    {

    }
}