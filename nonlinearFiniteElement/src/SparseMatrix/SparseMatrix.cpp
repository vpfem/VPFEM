#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(int& _N, int& _nnz, int* _row, int* _col, double* _value, bool& _type)
: m_sparsematrix{_N, _nnz, _row, _col, _value, _type}
{
}

SparseMatrix::~SparseMatrix()
{
}

// setters and getters
sp SparseMatrix::get_sp() const {return m_sparsematrix;};

void SparseMatrix::csrHost(int NumberCores)
{
    Convert2CSR(m_sparsematrix, 1);
}

void SparseMatrix::csrDevice()
{
    Convert2CSR(m_sparsematrix, 0);
}

// -- override the shift operator << 
std::ostream& operator<< (std::ostream &out, sp const& obj) 
{

    int iCounter = 0;
    for (int c = 0; c < obj.nnz; c++) {
        out << '\t' ;
        if (!obj.type) { // if the format is COO
        out << obj.row[c];
        } else if (c == obj.row[iCounter+1]-1) { // if the format is CSC
        out << obj.row[iCounter+1];
        iCounter++;
        }
        out << '\t' << obj.col[c] << ':' << '\t' << obj.value[c] << '\n';
    }
    return out;
};