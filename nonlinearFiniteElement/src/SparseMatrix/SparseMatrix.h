#pragma once 
#include <iostream>
#include <algorithm>
#include <cuda_runtime.h>

#include "sp.h"
#include "Convert2CSR.h"

class SparseMatrix {
public:

private:
    sp m_sparsematrix; // is a sparse matrix
    int* indices; // indices vector used for sorting and partitioning

public:
    SparseMatrix(int&, int&, int*, int*, double*, bool&);
    ~SparseMatrix();
    sp get_sp() const;
    void csrHost(int);
    void csrDevice();
private:

};
