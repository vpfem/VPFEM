#include "StiffnessMatrixGPU.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

StiffnessMatrixGPU::StiffnessMatrixGPU(double* mat, Geometry &geo, int n)
  : StiffnessMatrixFirstOrder(mat, geo, n)
{
  int device = -1;
  cudaGetDevice(&device);
  // copy from the material matarix
  //cudaMallocManaged(&D_d, 6*sizeof(double));
  //cudaMemcpy(D_d, material->materialMatrix, 6*sizeof(double), cudaMemcpyHostToDevice);
  cudaDeviceSynchronize();
  INFO("StiffnessMatrixGPU created by CPU");
};

StiffnessMatrixGPU::~StiffnessMatrixGPU()
{
  INFO("StiffnessMatrixGPU deleted by CPU");
  //cudaFree(D_d);
}

__global__ void constantCreatorKernel(int n, double* c, double* x, double* y, int* mesh, StiffnessMatrixGPU *s)
{
  //printf("in the function\n blockDim.x = %d, gridDim.x = %d, blockIdx.x = %d\n", blockDim.x,gridDim.x, blockIdx.x);
  int index = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride)
    {
      //printf("i is %d stride is %d threadID = %d\n",i,stride,threadIdx.x);
      s->constantCreator(i, c, x, y, mesh);
    }
};


__global__ void StiffnessMatrixKernel(  int n, int nip, double* in, int* ip, double* iw, double* c, double* D, int* mesh, double* k, int* i_index, int *j_index, int* dofFree, const double* thickness, StiffnessMatrixGPU *obj)
{
  int index  = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = index; i < n; i += stride)
    {
      obj->stiffnessMatrixCalculation(i, nip, in, ip, iw, c, D, mesh, k, i_index, j_index, dofFree, thickness);
    }
}

Sparse& StiffnessMatrixGPU::GetStiffnessMatrix()
{
  blockSize = 32;
  //numberOfElements=33;
  int numBlocks = (numberOfElements + blockSize-1)/blockSize;
  constantCreatorKernel<<<numBlocks, blockSize>>>(numberOfElements, c, geometry->get_x(), geometry->get_y(), geometry->get_mesh(), this);
  cudaDeviceSynchronize();
  numBlocks = (simulationSize + blockSize-1)/blockSize;
  //Timer timer("Time spend in GPU: ");
  StiffnessMatrixKernel<<<numBlocks, blockSize>>>(numberOfElements, nipSquared, integrationNode, integrationPos, integrationWeight, c, material, geometry->get_mesh(), stiffMat->value, stiffMat->i, stiffMat->j , geometry->get_Dof().get_free(), geometry->get_thickness(),this);
  //gpuErrchk( cudaPeekAtLastError() );
  //gpuErrchk( cudaDeviceSynchronize() );
  cudaDeviceSynchronize();
  return *stiffMat;
}

__global__ void DisplacementToStrainKernel(  int N, double * displacement, int * dofFree, int * mesh, double* strain,\
  int strainPerElement, int nip ,double* in, int* ip, double* c, StiffnessMatrixGPU *obj)
{
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i < N)
    obj->DisplacementToStrainEachElement(i, displacement, dofFree, mesh, strain, strainPerElement, nip, in, ip, c);
}
void StiffnessMatrixGPU::DisplacementToStrain(double* displacement, double* strain) 
{
  blockSize = 32;
  int numBlocks = (numberOfElements + blockSize-1)/blockSize;
  DisplacementToStrainKernel<<<numBlocks, blockSize>>>(numberOfElements, \
                            displacement, geometry->get_Dof().get_free(), \
                            geometry->get_mesh(), strain, 3, \
                            nipSquared, integrationNode, integrationPos, c, this);
  cudaDeviceSynchronize();
}