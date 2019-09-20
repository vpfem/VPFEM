#ifndef SP_H
#define SP_H

struct sp
{
    int N; // number of rows
    int nnz; // number of nonzeros
    int* row;
    int* col;
    double* value;
    bool type; // True -> CSR False -> COO
};

#endif