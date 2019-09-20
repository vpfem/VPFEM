#include "Convert2CSR.h"

void removeDuplicate(int i, int* row, int* col, double* value)
{
    int incr;
    int counter = row[i];
    while ( counter < row[i+1])
    {
        incr = 1;
        while (col[counter] == col[counter + incr] )
        {
            value[counter] += value[counter + incr];
            value[counter + incr] = 0;
            incr++;
        }
        counter += incr;
    }
}

__global__ void removeDuplicateKernel(int N, int* row, int* col, double* value)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N)
    {
        removeDuplicate(i, row, col, value);
    }
}

__global__ void deductOneColKernel(int nnz, int* col)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < nnz)
    {
        col[i]--;
    }
}

__global__ void reorderCopyKernel(int nnz, int* index, int* rowOld, int* colOld, double* valueOld, int* rowNew, int* colNew, double* valueNew)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < nnz)
    {
        //printf("i is %d and index is %d \n", i, index[i]);
        *(rowNew + *(index + i)) = *(rowOld + i);
        *(colNew + *(index + i)) = *(colOld + i);
        *(valueNew + *(index + i)) = *(valueOld + i);
    }
}

Convert2CSR::Convert2CSR(sp& A, int _device) 
: m_A(A), m_device(_device) 
{
    if (m_device == 0)
        Convert2CSRDevice();
    else
        Convert2CSRHost();
}

void Convert2CSR::Convert2CSRHost()
{
    cudaMallocManaged(&m_index, m_A.nnz*sizeof(int));
    cudaMallocManaged(&rowPtr, (m_A.N + 1)*sizeof(int));
    sortHost();
    reorderInplaceHost();
    csrRow();
    removeDuplicatesHost();
    partitionHost();
    reorderInplaceHost();
    csrRow();
    cudaFree(m_A.row);
    m_A.row = rowPtr;
    m_A.type = true;
    deductOneColHost();
}

void Convert2CSR::Convert2CSRDevice()
{
    cudaMallocManaged(&m_index, m_A.nnz*sizeof(int));
    cudaMallocManaged(&rowPtr, (m_A.N + 1)*sizeof(int));
    sortDevice();
    Timer timer("reorder time");
    timer.Start();
    //reorderInplaceDevice();
    reorderCopyDevice();
    timer.Stop();
    csrRow();
    removeDuplicatesDevice();
    partitionDevice();
    //reorderInplaceDevice();
    timer.Start();
    reorderCopyDevice();
    timer.Stop();
    csrRow();
    cudaFree(m_A.row);
    m_A.row = rowPtr;
    m_A.type = true;
    deductOneColDevice();
}

Convert2CSR::~Convert2CSR()
{
    cudaFree(m_index);
}

void Convert2CSR::partitionHost()
{
double* v = m_A.value;
for (int c = 0; c < m_A.nnz; c++) { m_index[c] = c;};
std::stable_partition(m_index, m_index + m_A.nnz, 
    [v](int i) 
    {
        return std::abs(v[i]) > 0.0001;
    });
}
    
void Convert2CSR::sortHost()
{
    for (int c = 0; c < m_A.nnz; c++) { m_index[c] = c;};
    std::sort(m_index, m_index + m_A.nnz, Compare(m_A.row,m_A.col));
}

void Convert2CSR::csrRow()
{
    int rowCounter = 0;
    for (int i = 0; i < m_A.N; i++)
    {
        while (*(m_A.row + rowCounter) == i+1)
        {
            rowCounter++;
        }
        rowPtr[i+1] = rowCounter;
    }
    rowPtr[0] = 0;
    m_A.nnz = rowPtr[m_A.N];
}
    
void Convert2CSR::removeDuplicatesHost()
{
    for (int i = 0; i < m_A.N; i++)
    {
        removeDuplicate(i, rowPtr, m_A.col, m_A.value);
    }
}

void Convert2CSR::reorderInplaceHost() 
{
    int index[m_A.nnz]; 
    for (int i=0; i< m_A.nnz; i++) {index[m_index[i]] = i;} // create array ordered indexes 
    for (int i=0; i< m_A.nnz; i++) {std::cout<<index[i]<<"\n";} // create array ordered indexes 

    for (int i=0; i < m_A.nnz; i++) 
    { 
        while (index[i] != i) 
        { 
            int  oldTargetI  = index[index[i]]; 
            int  oldTargetR  = m_A.row[index[i]]; 
            int  oldTargetC  = m_A.col[index[i]]; 
            double  oldTargetV  = m_A.value[index[i]]; 
            
            // Place row[i] at its target (or correct) 
            // position. Also copy corrected index for 
            // new position 
            m_A.row[index[i]] = m_A.row[i]; 
            m_A.col[index[i]] = m_A.col[i]; 
            m_A.value[index[i]] = m_A.value[i]; 
            index[index[i]] = index[i]; 
            
            // Copy old target values to row[i] and 
            // index[i] 
            index[i] = oldTargetI; 
            m_A.row[i]   = oldTargetR; 
            m_A.col[i]   = oldTargetC; 
            m_A.value[i]   = oldTargetV; 
        } 
    } 
}

void Convert2CSR::deductOneColHost()
{
    for (int i = 0; i < m_A.nnz; i++)
        m_A.col[i]--;
}

void Convert2CSR::partitionDevice()
{
    for (int c = 0; c < m_A.nnz; c++) { m_index[c] = c;};
    cudaDeviceSynchronize();
    thrust::stable_partition(thrust::device, m_index, m_index + m_A.nnz, Compare2zero(m_A.value));
}

void Convert2CSR::sortDevice()
{
    for (int c = 0; c < m_A.nnz; c++) { m_index[c] = c;};
    cudaDeviceSynchronize();
    thrust::sort(thrust::device, m_index, m_index + m_A.nnz, Compare(m_A.row,m_A.col));
}

void Convert2CSR::removeDuplicatesDevice()
{
    int blockSize = 32;
    int numBlocks = (m_A.N + blockSize-1)/blockSize;
    removeDuplicateKernel<<<numBlocks, blockSize>>>(m_A.N, rowPtr, m_A.col, m_A.value);
    cudaDeviceSynchronize();
}

void Convert2CSR::reorderInplaceDevice() 
{
    int* index; cudaMallocManaged(&index, m_A.nnz*sizeof(int));
    for (int i=0; i< m_A.nnz; i++) {index[m_index[i]] = i;} // create array ordered indexes 
    cudaDeviceSynchronize();
    thrust::sort_by_key(thrust::device, index, index + m_A.nnz, m_A.row);
    for (int i=0; i< m_A.nnz; i++) {index[m_index[i]] = i;} // create array ordered indexes 
    cudaDeviceSynchronize();
    thrust::sort_by_key(thrust::device, index, index + m_A.nnz, m_A.col);
    for (int i=0; i< m_A.nnz; i++) {index[m_index[i]] = i;} // create array ordered indexes 
    cudaDeviceSynchronize();
    thrust::sort_by_key(thrust::device, index, index + m_A.nnz, m_A.value);
    cudaFree(index);
}

void Convert2CSR::reorderCopyDevice()
{
    int blockSize = 32;
    int numBlocks = (m_A.nnz + blockSize-1)/blockSize;
    int* index; cudaMallocManaged(&index, m_A.nnz*sizeof(int)); for (int i=0; i< m_A.nnz; i++) {index[m_index[i]] = i;} // create array ordered indexes 
    int* t_row; cudaMallocManaged(&t_row, m_A.nnz*sizeof(int));
    int* t_col; cudaMallocManaged(&t_col, m_A.nnz*sizeof(int));
    double* t_val; cudaMallocManaged(&t_val, m_A.nnz*sizeof(double));
    reorderCopyKernel<<<numBlocks, blockSize>>>(m_A.nnz, index, m_A.row, m_A.col, m_A.value, t_row, t_col, t_val);
    cudaFree(m_A.row);
    cudaFree(m_A.col);
    cudaFree(m_A.value);
    m_A.row = t_row;
    m_A.col = t_col;
    m_A.value = t_val;
}

void Convert2CSR::deductOneColDevice()
{
    int blockSize = 32;
    int numBlocks = (m_A.nnz + blockSize-1)/blockSize;
    deductOneColKernel<<<numBlocks, blockSize>>>(m_A.nnz, m_A.col);
}

// -- Compare struct 
Compare::Compare(int* _row, int* _col)
  : row(_row), col(_col) 
{ };

bool Compare::operator()(int i, int j) 
{
    if (row[i] == 0 )
        return false;
    if (row[j] == 0)
        return true;
    if (row[i] == row[j])
        return col[i] < col[j];
    return row[i] < row[j];
}

// -- Compare struct 
Compare2zero::Compare2zero(double* value)
  : v(value)
{ };

bool Compare2zero::operator()(int i) 
{
    return abs(v[i]) > 0.0001;
}