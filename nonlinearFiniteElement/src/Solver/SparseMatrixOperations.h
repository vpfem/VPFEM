#pragma once

#include <iostream>
#include <algorithm>

template <typename INDEX, typename VALUETYPE>
class SparseMatrixOperations{
public:
    SparseMatrixOperations(){

    }

    // performs matrix vector multiplication on CSR matrix
    // matrix should be N*N
    // y = A*x
    void csrmv(int N, const INDEX* rowPtrA, const INDEX* colIndA, const VALUETYPE* valueA, const VALUETYPE* x, VALUETYPE* y)
    {
        const VALUETYPE* valueA_end = valueA + *(rowPtrA + N);
        const INDEX* colIndA_start = colIndA;
        std::fill(y, y+N, 0);
        while (valueA != valueA_end)
        {
            *y = *y + *(valueA) * (*(x + *colIndA));
            colIndA++;
            if (colIndA == colIndA_start + *(rowPtrA+1)) 
            { 
                rowPtrA++;
                y++;
            }
            valueA++;
        }
    };

    // copy
    // y = x
    void copy (int N, const VALUETYPE* x, VALUETYPE* y)
    {
        const VALUETYPE* x_end = x + N;
        while (x != x_end)
        {
            *y = *x;
            x++;
            y++;
        }
    }


    // This is equivalent to cublas<t>dot()
    // dot product of two vectors
    void dot(int N, const VALUETYPE* x, VALUETYPE* y, VALUETYPE* result)
    {
        *result = 0.0;
        auto x_end = x + N;
        while (x != x_end)
        {
            *result = *result + (*x)*(*y);
            x++;
            y++;
        }
    };

    // This is equivalent to cublas<t>scal()
    void scale(int N, const VALUETYPE alpha, VALUETYPE* x)
    {
        auto x_end = x + N;
        while (x != x_end)
        {
            *x = alpha*(*x);
            x++;
        }
    };

    // This is equivalent to cublas<t>axpy()
    // y = alpha * x
    void axpy(int N, const VALUETYPE alpha, const VALUETYPE* x, VALUETYPE* y)
    {
        auto x_end = x + N;
        while (x != x_end)
        {
            *y = alpha*(*x) + *y;
            x++;
            y++;
        }
    };
};
 