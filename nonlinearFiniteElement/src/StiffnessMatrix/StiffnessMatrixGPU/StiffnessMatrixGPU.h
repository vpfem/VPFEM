#ifndef STIFFNESSMATRIXGPU_H
#define STIFFNESSMATRIXGPU_H

#include <iostream>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "../StiffnessMatrixFirstOrder/StiffnessMatrixFirstOrder.h"
#include "../../Log/Log.h"
#include "../../Timer/Timer.h"
#include "../../Sparse/Sparse.h"

class StiffnessMatrixGPU : public StiffnessMatrixFirstOrder
{
  // variables
public:
  double* D_d; // required arrays to copy to GPU
private:
    int blockSize;
  // methods
public:
  StiffnessMatrixGPU(double*, Geometry&,   int);
  ~StiffnessMatrixGPU();
  Sparse& GetStiffnessMatrix() override;
  void DisplacementToStrain(double* displacement, double* strain);
private:

};
#endif
