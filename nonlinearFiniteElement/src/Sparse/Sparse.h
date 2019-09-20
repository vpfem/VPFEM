#ifndef SPARSE_H
#define SPARSE_H

#include <iostream>
#include <algorithm> 
#include <thrust/sort.h>

#include "../Log/Log.h"
#include "../Timer/Timer.h"
#include "SparseMatrix/SparseMatrix.h"
#include "SparseMatrix/Convert2CSR.h"
#include "SparseMatrix/sp.h"
//#include "Assembly.h"

//macro
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#endif

//#include <cuda_runtime_api.h>
//#include <cuda.h>
//#include  <cusolverDn.h>

class Sparse
{
  //variables
public:
  double* value;
  int* i;
  int* j;
  int valueSize; //how many nonzeros
  int type;
private:
  int numberOfRows; //number of rows
  int numberOfColumns; //number of columns
  bool symmetry = false; //by default matrix is not symmetry
  //cuda inverse related
  
  // Methods
public:
  Sparse(int, int, int);
  Sparse(int, int);
  Sparse();
  ~Sparse();
  void STLAssemble(double);
  void STLAssemble2(double);
  void STLAssemble3(double);
  void ThrustAssemble(double);
  void ThrustAssemble2(double);
  void ThrustAssemble3(double);
  
  void set_numberOfRows(int);
  void set_numberOfColumns(int);
  void set_valueSize(int);
  void set_i(int*);
  void set_j(int*);
  void set_x(double*);
  int get_valueSize() const;
  int get_numberOfRows() const;
  int get_numberOfColumns() const;
  double * get_value() const;
  int * get_i() const;
  int * get_j() const;
  int get_type() const;
  
private:
  void Print();

  // staic methods
public:
  static int Printer(Sparse & s);
  static double Compare(Sparse& , Sparse&);
};

std::ostream& operator<< (std::ostream &, Sparse const& );


struct sort_indices {
private:
  int* dofSorted; // the variable that index is going to sorted base of
public:
  CUDA_HOSTDEV sort_indices(int*);
  CUDA_HOSTDEV bool operator()(int i, int j) const;
};

struct sort_indices_j {
private:
  double* x; // holds values of stiffness matrix
  int* dofSorted; // the variable that index is going to sorted base of
public:
  CUDA_HOSTDEV sort_indices_j(int*, double*); 
  CUDA_HOSTDEV bool operator()(int i, int j);
};


// -- new Sort
struct newSort {
  int nRow;
  int* i; int* j;
  unsigned long long int a;
  unsigned long long int b;
  CUDA_HOSTDEV newSort(int*, int*, int);
  CUDA_HOSTDEV bool operator()(int, int);
};
#endif
