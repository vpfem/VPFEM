#include "Sparse.h"

Sparse::Sparse() {};

Sparse::Sparse(int x_size, int rowSize, int columnSize)
  : valueSize(x_size), numberOfRows(rowSize), numberOfColumns(columnSize)
{
  INFO("Sparse created");
  cudaMallocManaged(&i, valueSize*sizeof(int));
  cudaMallocManaged(&j, valueSize*sizeof(int));  
  cudaMallocManaged(&value, valueSize*sizeof(double));
  cudaMemset(j,0,valueSize*sizeof(int));
  cudaMemset(i,0,valueSize*sizeof(int));
  cudaMemset(value, 0 , valueSize*sizeof(double));
  type = 0; // means the sparse matrix is COO format
}

Sparse::Sparse(int x_size, int sizeOfMatrix)
  : Sparse(x_size, sizeOfMatrix, sizeOfMatrix)
{
  symmetry = true;
};

Sparse::~Sparse() {
  INFO("Sparse deleted");
  cudaFree(i);
  cudaFree(j);
  cudaFree(value);
}

void Sparse::STLAssemble3(double tol)
{
  sp matrix{numberOfRows, valueSize, i, j, value, type};
  Convert2CSR(matrix, 1);
  valueSize = matrix.nnz;
  type = matrix.type;
  i = matrix.row;
  j = matrix.col;
  
}
void Sparse::ThrustAssemble3(double tol)
{
  sp matrix{numberOfRows, valueSize, i, j, value, type};
  Convert2CSR(matrix, 0);
  valueSize = matrix.nnz;
  type = matrix.type;
  i = matrix.row;
  j = matrix.col;
  value = matrix.value;
}

// This function return a CSC matrix
void Sparse::STLAssemble(double tol)
// coo -> a triplet sparse format
{
  type = 1; // type is CSC
  int * a_temp_i;  int * a_temp_j; double * a_temp_value;
  // build the indices vector
  int* indices = new int[valueSize];
  for ( int c = 0; c < valueSize; c++) { indices[c] = c;};
  //// timer
    //Timer *timer = new Timer("Time spend in sorting for STL Assembley: ");
  // sorting i
  std::sort(indices,indices+valueSize, sort_indices(i));
  // number on non-zeros in each row
  int* nnz_inRow = new int[numberOfRows+1](); 
  int rowCounter = 0;
  for (int c = 0; c < valueSize; c++) {
    if (i[indices[c]] == rowCounter) {
      nnz_inRow[rowCounter] = c+1;
    } else {
      rowCounter++;
    }
  }
  // sorting j 
  cudaMallocManaged(&a_temp_i, (numberOfRows+1)*sizeof(int));
  a_temp_i[0] = 0;// rowPtr
  for (int c = 1; c <= numberOfRows; c++) {
    std::sort(indices + nnz_inRow[c-1], indices + nnz_inRow[c], sort_indices_j(j,value));
    // -- remove the zero entries (counting the number of nnz in each row)
    a_temp_i[c] = 0;
    for (int count = nnz_inRow[c-1]; count < nnz_inRow[c]; count++) 
      if (abs(value[indices[count]]) > tol) {a_temp_i[c]++;}
    a_temp_i[c] += a_temp_i[c-1];
  }
  //delete timer; 
  //Timer t; t.Start();
  // copy to new array
  cudaMallocManaged(&a_temp_j, (a_temp_i[numberOfRows])*sizeof(int)); 
  cudaMallocManaged(&a_temp_value, (a_temp_i[numberOfRows])*sizeof(double));
  int vs = 0; // new value size
  for (int c = 0; c < valueSize; c++) {
    if (abs(value[indices[c]]) > tol) {
      //a_temp_i[vs] = i[indices[c]];
      a_temp_j[vs] = j[indices[c]] - 1;
      a_temp_value[vs] = value[indices[c]];
      vs++;
    }
  }
  valueSize = a_temp_i[numberOfRows];
  cudaFree(i); cudaFree(j); cudaFree(value);
  i = a_temp_i;
  j = a_temp_j;
  value = a_temp_value;  
  delete[] indices;
  delete[] nnz_inRow;
  //t.Stop();
}

// This function return a CSC matrix using thrust algorithm
void Sparse::ThrustAssemble(double tol)
// coo -> a triplet sparse format
{
  cudaDeviceSynchronize();
  type = 1; // type is CSC
  int * a_temp_i;  int * a_temp_j; double * a_temp_value;
  // build the indices vector
  int* indices;  cudaMallocManaged(&indices, valueSize*sizeof(int));
  for ( int c = 0; c < valueSize; c++) { indices[c] = c;};
  //// timer
    //Timer *timer = new Timer("Time spend in sorting for Thrust Assembley: ");
  // sorting i
  thrust::sort(indices,indices+valueSize, sort_indices(i));
  // number on non-zeros in each row
  int* nnz_inRow;
  cudaMallocManaged(&nnz_inRow, (numberOfRows+1)*sizeof(int));
  cudaMemset(nnz_inRow, 0, (numberOfRows+1)*sizeof(int));
  int rowCounter = 0;
  for (int c = 0; c < valueSize; c++) {
    if (i[indices[c]] == rowCounter) {
      nnz_inRow[rowCounter] = c+1;
    } else {
      rowCounter++;
    }
  }
  // sorting j 
  cudaMallocManaged(&a_temp_i, (numberOfRows+1)*sizeof(int));
  a_temp_i[0] = 0;// rowPtr
  for (int c = 1; c <= numberOfRows; c++) {
    thrust::sort(indices + nnz_inRow[c-1], indices + nnz_inRow[c], sort_indices_j(j,value));
    // -- remove the zero entries (counting the number of nnz in each row)
    a_temp_i[c] = 0;
    for (int count = nnz_inRow[c-1]; count < nnz_inRow[c]; count++) 
      if (abs(value[indices[count]]) > tol) {a_temp_i[c]++;}
    a_temp_i[c] += a_temp_i[c-1];
  }
  //delete timer; Timer t("Time for rest of  Assembley: ");
  // copy to new array
  cudaMallocManaged(&a_temp_j, (a_temp_i[numberOfRows])*sizeof(int)); 
  cudaMallocManaged(&a_temp_value, (a_temp_i[numberOfRows])*sizeof(double));
  int vs = 0; // new value size
  for (int c = 0; c < valueSize; c++) {
    if (abs(value[indices[c]]) > tol) {
      //a_temp_i[vs] = i[indices[c]];
      a_temp_j[vs] = j[indices[c]] - 1;
      a_temp_value[vs] = value[indices[c]];
      vs++;
    }
  }
  valueSize = a_temp_i[numberOfRows];
  cudaFree(i); cudaFree(j); cudaFree(value);
  i = a_temp_i;
  j = a_temp_j;
  value = a_temp_value;  
  cudaFree(indices);
  cudaFree(nnz_inRow);
}

double Sparse::Compare(Sparse& A, Sparse& B)
// return the max Error between the two matrices
{
  double MaxError = 0.00;
  for (int counter = 0; counter < A.valueSize; counter++)
    MaxError = fmax(MaxError, fabs(A.value[counter] - B.value[counter]));
  return MaxError;
}
// -- new assemblers
void Sparse::STLAssemble2(double tol) {
  
  type = 1; // type is CSC
  // build the indices vector
  int* indices;  cudaMallocManaged(&indices, valueSize*sizeof(int));
  for ( int c = 0; c < valueSize; c++) { indices[c] = c;};
  //// timer
  //Timer *timer = new Timer("Time spend in sorting for STL Assembley2: ");
  // sorting i
  std::sort(indices,indices+valueSize, newSort(i,j,numberOfRows));
  //thrust::sort(thrust::device,indices,indices+valueSize, newSort(i,j,numberOfRows));
  //delete timer; timer = new Timer("Time for rest of STL Assembley2: ");
  // -- define rowPtr
  int* a_temp_i; cudaMallocManaged(&a_temp_i, (numberOfRows+1)*sizeof(int));
  cudaMemset(a_temp_i,0,(numberOfRows+1)*sizeof(int));
  // -- find how many zeros are present
  size_t zeroCounter = 0;
  while (!i[indices[zeroCounter]]) {zeroCounter++;}
  // -- add duplicates and creat the rowPtr
  size_t count;
  for (int c = 0; c < valueSize-1; c++) {
    count = c+1;
    while (i[indices[c]] == i[indices[count]] && j[indices[c]] == j[indices[count]]) {
      value[indices[c]] += value[indices[count]];
      i[indices[count]] = 0; j[indices[count]] = 0; value[indices[count]] = 0;
	count++;
    }
    if (abs(value[indices[c]]) > tol) {
      a_temp_i[i[indices[c]]]++;
    } 
    c = count-1;
  }
  if (abs(value[indices[valueSize-1]]) > tol) {a_temp_i[i[indices[valueSize-1]]]++;}
  // -- print to screen
  //for (int c = 0; c < valueSize; c++)
  //  std::cout<< c << '\t' << i[indices[c]] << '\t' << j[indices[c]] << '\t' << value[indices[c]] << '\n';
  // -- sum rowPtr
  for (int c = 1; c<=numberOfRows; c++)
    a_temp_i[c] += a_temp_i[c-1];
  //delete timer; Timer t("Time for copy of STL Assembley2: ");
  // copy to new variables
  int* a_temp_j; cudaMallocManaged(&a_temp_j,     (a_temp_i[numberOfRows])*sizeof(int)); 
  double*   a_temp_value; cudaMallocManaged(&a_temp_value, (a_temp_i[numberOfRows])*sizeof(double));
  int valueCounter = 0;
  for (int c = 0; c < valueSize; c++) {
    if(abs(value[indices[c]]) > tol) {
      a_temp_j[valueCounter] = j[indices[c]]-1;
      a_temp_value[valueCounter] = value[indices[c]];
      valueCounter++;
    }
  }
  // change variables
  valueSize = a_temp_i[numberOfRows];
  cudaFree(i); cudaFree(j); cudaFree(value);
  i = a_temp_i;
  j = a_temp_j;
  value = a_temp_value;  
  cudaFree(indices);
}

// -- new assemblers thrust
void Sparse::ThrustAssemble2(double tol) {
  type = 1; // type is CSC
  // build the indices vector
  int* indices;  cudaMallocManaged(&indices, valueSize*sizeof(int));
  for ( int c = 0; c < valueSize; c++) { indices[c] = c;};
  //// timer
  //Timer *timer = new Timer("Time spend in sorting for STL Assembley2: ");
  // sorting i
  //std::sort(indices,indices+valueSize, newSort(i,j,numberOfRows));
  cudaDeviceSynchronize();
  thrust::sort(thrust::device,indices,indices+valueSize, newSort(i,j,numberOfRows));
  //delete timer; Timer t("Time for rest of STL Assembley2: ");
  // -- define rowPtr
  int* a_temp_i; cudaMallocManaged(&a_temp_i, (numberOfRows+1)*sizeof(int));
  cudaMemset(a_temp_i,0,(numberOfRows+1)*sizeof(int));
  // -- find how many zeros are present
  size_t zeroCounter = 0;
  while (!i[indices[zeroCounter]]) {zeroCounter++;}
  // -- add duplicates and creat the rowPtr
  size_t count;
  for (int c = zeroCounter; c < valueSize-1; c++) {
    count = c+1;
    while (i[indices[c]] == i[indices[count]] && j[indices[c]] == j[indices[count]]) {
      value[indices[c]] += value[indices[count]];
      i[indices[count]] = 0; j[indices[count]] = 0; value[indices[count]] = 0;
	count++;
    }
    if (abs(value[indices[c]]) > tol) {
      a_temp_i[i[indices[c]]]++;
    } else {
      i[indices[c]] = 0; 
    }
    c = count-1;
  }
  if (abs(value[indices[valueSize-1]]) > tol) {a_temp_i[i[indices[valueSize-1]]]++;}
  // -- sum rowPtr
  for (int c = 1; c<=numberOfRows; c++)
    a_temp_i[c] += a_temp_i[c-1];
  // copy to new variables
  int* a_temp_j; cudaMallocManaged(&a_temp_j,     (a_temp_i[numberOfRows])*sizeof(int)); 
  double*   a_temp_value; cudaMallocManaged(&a_temp_value, (a_temp_i[numberOfRows])*sizeof(double));
  size_t valueCounter = 0;
  for (int c = zeroCounter; c < valueSize; c++) {
    if(i[indices[c]]) {
      a_temp_j[valueCounter] = j[indices[c]]-1;
      a_temp_value[valueCounter] = value[indices[c]];
      valueCounter++;
    }
  }
  // change variables
  valueSize = a_temp_i[numberOfRows];
  cudaFree(i); cudaFree(j); cudaFree(value);
  i = a_temp_i;
  j = a_temp_j;
  value = a_temp_value;  
  cudaFree(indices);
}

// Setters
void Sparse::set_i(int* index_i) { i = index_i;}
void Sparse::set_j(int* index_j) { j = index_j;}
void Sparse::set_x(double* x) { value = x;}
void Sparse::set_numberOfRows(int x) { numberOfRows = x;}
void Sparse::set_valueSize(int x) { valueSize = x;}
void Sparse::set_numberOfColumns(int x) {
  numberOfColumns = x;
  if (x == 1 ) { // if we are dealing with vecotrs 
    j = new int[valueSize];
    for (int ci = 0; ci < valueSize; ci++) {
      j[ci] = 1;
    }
  }
}
  
// Geters
int Sparse::get_valueSize() const {return valueSize;};
int Sparse::get_numberOfRows() const { return numberOfRows;}
int Sparse::get_numberOfColumns() const {return numberOfColumns;}
double * Sparse::get_value() const {return value;}
int * Sparse::get_i() const {return i;}
int * Sparse::get_j() const {return j;}
int Sparse::get_type() const {return type;}

// -- printers
void Sparse::Print() {
  std::cout << "\033[1;34m[SparseMatrix]: \033[0m" << numberOfRows  << " x " << numberOfColumns << " ";
  if (symmetry)
    std::cout << "\033[32msymmetry\033[0m " ;
  int size = valueSize <72 ? valueSize : 72; 
  std::cout << "print size: " << size << ", nnz = " << valueSize << std::endl;
  for (int counter = 0; counter < size; counter++)
    std::cout << i[counter] << "\t"<< j[counter] <<"\t: " << value[counter] <<std::endl;
};

int Sparse::Printer(Sparse& s) {
  s.Print();
  return 0;
};

// -- override the cout << oprator 
std::ostream& operator<< (std::ostream &out, Sparse const& sp) {
  const double* x = sp.get_value();
  const int* j = sp.get_j();
  const int* i = sp.get_i();
  int iCounter = 0;
  for (int c = 0; c < sp.get_valueSize(); c++) {
    out << '\t' ;
    if (!sp.get_type()) { // if the format is COO
      out << i[c];
	} else if (c == i[iCounter+1]-1) { // if the format is CSC
      out << i[iCounter+1];
      iCounter++;
    }
    out << '\t' << j[c] << ':' << '\t' << x[c] << '\n';
  }
  return out;
}

// ---------------------- sort_indices struct -------------------------
sort_indices::sort_indices(int* var)
  : dofSorted(var) {};

bool sort_indices::operator()(int i, int j) const { return dofSorted[i] < dofSorted[j];};

// ---------------------- sort_indices_j struct ------------------------- 

sort_indices_j::sort_indices_j(int* j, double* value)
  : dofSorted(j), x(value) { };

bool sort_indices_j::operator()(int i, int j) {
  if (dofSorted[i] == dofSorted[j]) {
    x[i] = x[i] + x[j]; x[j] = 0;
    return false;
  } else {return dofSorted[i] < dofSorted[j];}
};

// ----------------------- new Sort struct ---------------------
newSort::newSort(int* i_index, int* j_index, int numberOfRows)
  : i(i_index), j(j_index), nRow(numberOfRows) {
};

bool newSort::operator()(int first, int second) {
  // compares the actual place in big stiffness matrix
  a = (unsigned long long int)(i[first]) *nRow + (unsigned long long int)(j[first]);
  b = (unsigned long long int)(i[second])*nRow + (unsigned long long int)(j[second]);
  return a < b;
}
