#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>

#include "../Log/Log.h"
#include "../Timer/Timer.h"
#include "../Sparse/Sparse.h"

//macro
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

class Assembly
{
public:
    Sparse &stiffMat;
    double tolorance;

    thrust::device_vector<double> d_value;
    thrust::device_vector<size_t> d_dofi;
    thrust::device_vector<size_t> d_dofj;

    
    thrust::device_vector<size_t> d_indices;

public:
    Assembly(Sparse& s);
    virtual ~Assembly();

    virtual void sort() = 0;
    virtual void calculateAssembly() = 0;
};

std::ostream& operator<< (std::ostream &, Assembly const& );

#endif