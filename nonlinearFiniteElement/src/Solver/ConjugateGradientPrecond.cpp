#include "ConjugateGradientPrecond.h"
namespace GPU {
    ConjugateGradientPrecondLU::ConjugateGradientPrecondLU(int N, int nz, double * val, int *rowPtr, int *colIndex, double * x, double* rhs, const double tol)
    {
        /* Conjugate Gradient variable */
        const int max_iter = 10000;
        int k = 0;
        double alpha, beta; 
        double r1;
        double numerator, denominator, nalpha;

        alpha = 1.0;
        beta = 0.0;

        double *p, *y, *omega;
        cudaMallocManaged(&p, N*sizeof(double));
        cudaMallocManaged(&y, N*sizeof(double));
        cudaMallocManaged(&omega, N*sizeof(double));

        /* Conjugate Gradient Precond variable */
        int nzILU0 = 2*N-1;
        double *valsILU0; cudaMallocManaged(&valsILU0, nz*sizeof(double));
        double *zm1; cudaMallocManaged(&zm1, N*sizeof(double));
        double *zm2; cudaMallocManaged(&zm2, N*sizeof(double));
        double *rm2; cudaMallocManaged(&rm2, N*sizeof(double));
        const double floatone = 1.0;
        const double floatzero = 0.0;
        /* Create CUBLAS context */
        cublasHandle_t cublasHandle = 0;
        cublasStatus_t cublasStatus;
        cublasStatus = cublasCreate(&cublasHandle);
        checkCudaErrors(cublasStatus);
        /* Create CUSPARSE context */
        cusparseHandle_t cusparseHandle = 0;
        cusparseStatus_t cusparseStatus;
        cusparseStatus = cusparseCreate(&cusparseHandle);
        
        /* create the analysis info object for the A matrix */
        cusparseSolveAnalysisInfo_t infoA = 0;
        cusparseStatus = cusparseCreateSolveAnalysisInfo(&infoA);

        /* Description of the A matrix*/
        cusparseMatDescr_t descr = 0;
        cusparseStatus = cusparseCreateMatDescr(&descr);

        /* Define the properties of the matrix */
        cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

        cusparseStatus = cusparseDcsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                N, nz, descr, val, rowPtr, colIndex, infoA);
        checkCudaErrors(cusparseStatus);
        /* Copy A data to ILU0 vals as input*/
        cudaMemcpy(valsILU0, val, nz*sizeof(double), cudaMemcpyDeviceToDevice);

        /* generate the Incomplete LU factor H for the matrix A using cudsparseDcsrilu0 */
        cusparseStatus = cusparseDcsrilu0(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, descr, valsILU0, rowPtr, colIndex, infoA);
        checkCudaErrors(cusparseStatus);
        /* Create info objects for the ILU0 preconditioner */
        cusparseSolveAnalysisInfo_t info_u;
        cusparseCreateSolveAnalysisInfo(&info_u);

        cusparseMatDescr_t descrL = 0;
        cusparseStatus = cusparseCreateMatDescr(&descrL);
        cusparseSetMatType(descrL,CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descrL,CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER);
        cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT); 

        cusparseMatDescr_t descrU = 0;
        cusparseStatus = cusparseCreateMatDescr(&descrU);
        cusparseSetMatType(descrL,CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descrU,CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER);
        cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT);
        cusparseDcsrsv_analysis(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, nz, descrU, val, rowPtr, colIndex, info_u);

        cublasDdot(cublasHandle, N, rhs, 1, rhs, 1, &r1); 

        while (r1 > tol*tol && k <= max_iter)
        {
            // Forward Solve, we can re-use infoA since the sparsity pattern of A matches that of L
            cusparseStatus = cusparseDcsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &floatone, descrL,
                                valsILU0, rowPtr, colIndex, infoA, rhs, y);
            checkCudaErrors(cusparseStatus);
            cusparseStatus = cusparseDcsrsv_solve(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, &floatone, descrU,
                                valsILU0, rowPtr, colIndex, info_u, y, zm1);
            checkCudaErrors(cusparseStatus);
            k++;

            if (k == 1)
            {
                cublasDcopy(cublasHandle, N, zm1, 1, p, 1);
            }
            else
            {
                cublasDdot(cublasHandle, N, rhs, 1, zm1, 1, &numerator);
                cublasDdot(cublasHandle, N, rm2, 1, zm2, 1, &denominator);
                beta = numerator/denominator;
                cublasDscal(cublasHandle, N, &beta, p, 1);
                cublasDaxpy(cublasHandle, N, &floatone, zm1, 1, p, 1) ;
            }

            cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nzILU0, &floatone, descrU, val, rowPtr, colIndex, p, &floatzero, omega);
            cublasDdot(cublasHandle, N, rhs, 1, zm1, 1, &numerator);
            cublasDdot(cublasHandle, N, p, 1, omega, 1, &denominator);
            alpha = numerator / denominator;
            cublasDaxpy(cublasHandle, N, &alpha, p, 1, x, 1);
            cublasDcopy(cublasHandle, N, rhs, 1, rm2, 1);
            cublasDcopy(cublasHandle, N, zm1, 1, zm2, 1);
            nalpha = -alpha;
            cublasDaxpy(cublasHandle, N, &nalpha, omega, 1, rhs, 1);
            cublasDdot(cublasHandle, N, rhs, 1, rhs, 1, &r1);
        }

        std::string report = "  iteration = " + std::to_string(k) + " residual = " + std::to_string(sqrt(r1));
        INFO(report);

        cusparseDestroySolveAnalysisInfo(infoA);
        cusparseDestroySolveAnalysisInfo(info_u);

        cusparseDestroy(cusparseHandle);
        cublasDestroy(cublasHandle);
        cudaFree(p);
        cudaFree(y);
        cudaFree(p);
        cudaFree(omega);
        cudaFree(valsILU0);
        cudaFree(zm1);
        cudaFree(zm2);
        cudaFree(rm2);
    }

    ConjugateGradientPrecondLU::~ConjugateGradientPrecondLU() { };

   ConjugateGradientPrecondChol::ConjugateGradientPrecondChol(int N, int nnz, double * val, int *rowPtr, int *colIndex, double * x, double* rhs, const double tol) 
    {
        /* Conjugate Gradient variable */
        const int max_iter = 10000;
        int k = 0;
        int m = N;
        double r1;
        double numerator, denominator, nalpha;
        double alpha = 1.0, beta = 0.0; 

        double *p, *y, *omega;
        cudaMallocManaged(&p, N*sizeof(double));
        cudaMallocManaged(&y, N*sizeof(double));
        cudaMallocManaged(&omega, N*sizeof(double));
        /* Conjugate Gradient Precond variable */
        double *valSV; cudaMallocManaged(&valSV, nnz*sizeof(double));
        double *zm1; cudaMallocManaged(&zm1, N*sizeof(double));
        double *zm2; cudaMallocManaged(&zm2, N*sizeof(double));
        double *rm2; cudaMallocManaged(&rm2, N*sizeof(double));
        const double floatone = 1.0;
        const double floatzero = 0.0;
        /* Create CUBLAS context */
        cublasHandle_t cublasHandle = 0;
        cublasStatus_t cublasStatus;
        cublasStatus = cublasCreate(&cublasHandle);
        checkCudaErrors(cublasStatus);
        /* Create CUSPARSE context */
        cusparseHandle_t cusparseHandle = 0;
        cusparseStatus_t cusparseStatus;
        cusparseStatus = cusparseCreate(&cusparseHandle);
        
        cusparseMatDescr_t descr_M = 0;
        cusparseMatDescr_t descr_L = 0;
        csric02Info_t info_M  = 0;
        csrsv2Info_t  info_L  = 0;
        csrsv2Info_t  info_Lt = 0;
        int pBufferSize_M;
        int pBufferSize_L;
        int pBufferSize_Lt;
        int pBufferSize;
        void *pBuffer = 0;
        int structural_zero;
        const cusparseSolvePolicy_t policy_M  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
        const cusparseSolvePolicy_t policy_L  = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
        const cusparseSolvePolicy_t policy_Lt = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
        const cusparseOperation_t trans_L  = CUSPARSE_OPERATION_NON_TRANSPOSE;
        const cusparseOperation_t trans_Lt = CUSPARSE_OPERATION_TRANSPOSE;

        // step 1: create a descriptor which contains
        // - matrix M is base-1
        // - matrix L is base-1
        // - matrix L is lower triangular
        // - matrix L has non-unit diagonal
        cusparseCreateMatDescr(&descr_M);
        cusparseSetMatIndexBase(descr_M, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatType(descr_M, CUSPARSE_MATRIX_TYPE_GENERAL);

        cusparseCreateMatDescr(&descr_L);
        cusparseSetMatIndexBase(descr_L, CUSPARSE_INDEX_BASE_ZERO);
        cusparseSetMatType(descr_L, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatFillMode(descr_L, CUSPARSE_FILL_MODE_LOWER);
        cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_NON_UNIT);

        // step 2: create a empty info structure
        // we need one info for csric02 and two info's for csrsv2
        cusparseCreateCsric02Info(&info_M);
        cusparseCreateCsrsv2Info(&info_L);
        cusparseCreateCsrsv2Info(&info_Lt);

        // step 3: query how much memory used in csric02 and csrsv2, and allocate the buffer
        cusparseDcsric02_bufferSize(cusparseHandle, m, nnz,
            descr_M, val, rowPtr, colIndex, info_M, &pBufferSize_M);
        cusparseDcsrsv2_bufferSize(cusparseHandle, trans_L, m, nnz,
            descr_L, val, rowPtr, colIndex, info_L, &pBufferSize_L);
        cusparseDcsrsv2_bufferSize(cusparseHandle, trans_Lt, m, nnz,
            descr_L, val, rowPtr, colIndex, info_Lt,&pBufferSize_Lt);

        pBufferSize = std::max(pBufferSize_M, std::max(pBufferSize_L, pBufferSize_Lt));

        // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
        cudaMalloc((void**)&pBuffer, pBufferSize);

        // step 4: perform analysis of incomplete Cholesky on M
        //         perform analysis of triangular solve on L
        //         perform analysis of triangular solve on L'
        // The lower triangular part of M has the same sparsity pattern as L, so
        // we can do analysis of csric02 and csrsv2 simultaneously.
        cusparseDcsric02_analysis(cusparseHandle, m, nnz, descr_M,
                                    val, rowPtr, colIndex, info_M,
                                    policy_M, pBuffer);
        cusparseStatus = cusparseXcsric02_zeroPivot(cusparseHandle, info_M, &structural_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == cusparseStatus){
        printf("A(%d,%d) is missing\n", structural_zero, structural_zero);
        }

        cusparseDcsrsv2_analysis(cusparseHandle, trans_L, m, nnz, descr_L,
            val, rowPtr, colIndex,
            info_L, policy_L, pBuffer);

        cusparseDcsrsv2_analysis(cusparseHandle, trans_Lt, m, nnz, descr_L,
            val, rowPtr, colIndex,
            info_Lt, policy_Lt, pBuffer);

        // step 5: M = L * L'
        cudaMemcpy(valSV, val, nnz*sizeof(double), cudaMemcpyDeviceToDevice);
        cusparseDcsric02(cusparseHandle, m, nnz, descr_M,
            valSV, rowPtr, colIndex, info_M, policy_M, pBuffer);
        /*
        cusparseStatus = cusparseXcsric02_zeroPivot(cusparseHandle, info_M, &numerical_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == cusparseStatus){
        printf("L(%d,%d) is zero\n", numerical_zero, numerical_zero);
        }
        */
        cublasDdot(cublasHandle, N, rhs, 1, rhs, 1, &r1);

        while (r1 > tol*tol && k <= max_iter)
        {
            // Forward Solve, we can re-use infoA since the sparsity pattern of A matches that of L
            cusparseStatus = cusparseDcsrsv2_solve(cusparseHandle, trans_L, m, nnz, &floatone, descr_L,
                                valSV, rowPtr, colIndex, info_L,
                                rhs, y, policy_L, pBuffer);
            checkCudaErrors(cusparseStatus);
            cusparseStatus = cusparseDcsrsv2_solve(cusparseHandle, trans_Lt, m, nnz, &floatone, descr_L,
                                valSV, rowPtr, colIndex, info_Lt,
                                y, zm1, policy_Lt, pBuffer);
            checkCudaErrors(cusparseStatus);
            k++;

            if (k == 1)
            {
                cublasDcopy(cublasHandle, N, zm1, 1, p, 1);
            }
            else
            {
                cublasDdot(cublasHandle, N, rhs, 1, zm1, 1, &numerator);
                cublasDdot(cublasHandle, N, rm2, 1, zm2, 1, &denominator);
                beta = numerator/denominator;
                cublasDscal(cublasHandle, N, &beta, p, 1);
                cublasDaxpy(cublasHandle, N, &floatone, zm1, 1, p, 1) ;
            }

            cusparseDcsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nnz, &floatone, descr_L, val, rowPtr, colIndex, p, &floatzero, omega);
            cublasDdot(cublasHandle, N, rhs, 1, zm1, 1, &numerator);
            cublasDdot(cublasHandle, N, p, 1, omega, 1, &denominator);
            alpha = numerator / denominator;
            cublasDaxpy(cublasHandle, N, &alpha, p, 1, x, 1);
            cublasDcopy(cublasHandle, N, rhs, 1, rm2, 1);
            cublasDcopy(cublasHandle, N, zm1, 1, zm2, 1);
            nalpha = -alpha;
            cublasDaxpy(cublasHandle, N, &nalpha, omega, 1, rhs, 1);
            cublasDdot(cublasHandle, N, rhs, 1, rhs, 1, &r1);
        }

        std::string report = "  iteration = " + std::to_string(k) + " residual = " + std::to_string(sqrt(r1));
        INFO(report);

        cudaFree(pBuffer);
        cusparseDestroyMatDescr(descr_M);
        cusparseDestroyMatDescr(descr_L);
        cusparseDestroyCsric02Info(info_M);
        cusparseDestroyCsrsv2Info(info_L);
        cusparseDestroyCsrsv2Info(info_Lt);

        cusparseDestroy(cusparseHandle);
        cublasDestroy(cublasHandle);
        cudaFree(p);
        cudaFree(y);
        cudaFree(p);
        cudaFree(omega);
        cudaFree(valSV);
        cudaFree(zm1);
        cudaFree(zm2);
        cudaFree(rm2);

    }
 
    ConjugateGradientPrecondChol::~ConjugateGradientPrecondChol() { }
}