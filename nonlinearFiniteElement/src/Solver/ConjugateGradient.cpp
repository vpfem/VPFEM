#include "ConjugateGradient.h"


namespace GPU {
    ConjugateGradient::ConjugateGradient(int N, int nz, double * val, int *rowPtr, int *colIndex, double * x, double* rhs, const double tol)
    {
        const int max_iter = 10000;
        int k = 0;
        double alpha, nalpha,beta, dot;
        double r0, r1;
        double *p, *omega;

        const double floatone = 1.0;
        const double floatzero = 0.0;


        alpha = 1.0;
        beta = 0.0;
        r0 = 0.;

        cudaMallocManaged(&p, N*sizeof(double));
        cudaMallocManaged(&omega, N*sizeof(double));

        /* Get handle to the CUBLAS context */
        cublasHandle_t cublasHandle = 0;
        cublasStatus_t cublasStatus;
        cublasStatus = cublasCreate(&cublasHandle);

        /* Get handle to the CUSPARSE context */
        cusparseHandle_t cusparseHandle = 0;
        cusparseStatus_t cusparseStatus;
        cusparseStatus = cusparseCreate(&cusparseHandle);

        cusparseMatDescr_t descr = 0;
        cusparseStatus = cusparseCreateMatDescr(&descr);

        cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

        cublasDdot(cublasHandle, N, rhs, 1, rhs, 1, &r1);
        while (r1 > tol*tol && k <= max_iter)
        {
            k++;

            if (k == 1)
            {
                cublasDcopy(cublasHandle, N, rhs, 1, p, 1);
            }
            else
            {
                beta = r1/r0;
                cublasDscal(cublasHandle, N, &beta, p, 1);
                cublasDaxpy(cublasHandle, N, &floatone, rhs, 1, p, 1) ;
            }
            cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &floatone, descr, val, rowPtr, colIndex, p, &floatzero, omega);
            cublasDdot(cublasHandle, N, p, 1, omega, 1, &dot);
            alpha = r1/dot;
            cublasDaxpy(cublasHandle, N, &alpha, p, 1, x, 1);
            nalpha = -alpha;
            cublasDaxpy(cublasHandle, N, &nalpha, omega, 1, rhs, 1);
            r0 = r1;
            cublasDdot(cublasHandle, N, rhs, 1, rhs, 1, &r1);
        }
        std::string report = "  iteration = " + std::to_string(k) + " residual = " + std::to_string(sqrt(r1));
        INFO(report);

        cusparseDestroy(cusparseHandle);
        cublasDestroy(cublasHandle);
        cudaFree(p);
        cudaFree(omega);
    }

    ConjugateGradient::~ConjugateGradient() {
    }
}

namespace SingleCPU {
    ConjugateGradient::ConjugateGradient(int N, int nz, double * val, int *rowPtr, int *colIndex, double * x, double* rhs, const double tol)
    {
        SparseMatrixOperations<int, double> sp;
        const int max_iter = 10000;
        int k = 0;
        double r0 = 0.0, r1, dot, alpha, nalpha;
        double beta = 0.0;
        double doubleOne = 1.0;
        double* p = new double[N];
        double* omega = new double[N];
        // dot prodoct
        sp.dot(N,rhs,rhs, &r1);
        while (r1 > tol*tol && k <= max_iter)
        {
            k++;
            if (k == 1) { sp.copy(N, rhs, p); }
            else 
            {
                beta = r1/r0;
                sp.scale(N,beta,p);
                sp.axpy(N, doubleOne, rhs, p);
            }
            sp.csrmv(N, rowPtr, colIndex, val, p, omega); 
            sp.dot(N,p,omega, &dot);
            alpha = r1/dot;
            sp.axpy(N, alpha, p, x);
            nalpha = -alpha;
            sp.axpy(N, nalpha, omega, rhs);
            r0 = r1;
            sp.dot( N, rhs, rhs, &r1);
        }

        std::string report = "  iteration = " + std::to_string(k) + " residual = " + std::to_string(std::sqrt(r1));
        INFO(report);

        // check the results
        //for (int i = 0; i < 4; i++) {std::cout << "F[" << i << "]= " << x[i] << std::endl;}
        // deletes
        delete[] p;
        delete[] omega;
    }

    ConjugateGradient::~ConjugateGradient() {
    }
}