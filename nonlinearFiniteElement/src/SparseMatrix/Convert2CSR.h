#ifndef Convert2CSR_H
#define Convert2CSR_H

#include <iostream>
#include <algorithm>
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

#include "Timer/Timer.h"
#include "sp.h"

//macro
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

class Convert2CSR 
{
private:
    sp& m_A;
    int *m_index;
    int *rowPtr;
    int m_device; // (0 - GPU) (1 - Single CPU) (2 - parallel CPU)
public:
    Convert2CSR(sp& A, int device); 
    ~Convert2CSR();
    void Convert2CSRHost();
    void Convert2CSRDevice();

private:
    void csrRow();
    void sortHost();
    void reorderHost();
    void removeDuplicatesHost();
    void partitionHost();
    void sortDevice();
    void reorderInplaceDevice();
    void reorderCopyDevice();
    void removeDuplicatesDevice();
    void partitionDevice();
    void reorderInplaceHost();
    void deductOneColDevice();
    void deductOneColHost();
};

// -- Compare struct 
struct Compare {
    int* row;
    int* col;
    CUDA_HOSTDEV Compare(int*, int*);
    CUDA_HOSTDEV bool operator()(int, int);
};
// -- Compare struct 
struct Compare2zero {
    double* v;
    CUDA_HOSTDEV Compare2zero(double*);
    CUDA_HOSTDEV bool operator()(int);
};

CUDA_HOSTDEV void removeDuplicate(int i, int* row, int* col, double* value);

// -- override the shift operator << 
std::ostream& operator<< (std::ostream &out, sp const& obj);

#endif