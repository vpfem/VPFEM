#ifndef STIFFNESSMATRIXSINGLECPU
#define STIFFNESSMATRIXSINGLECPU

#include <iostream>
#include <thread>
#include <string>
#include "../StiffnessMatrixFirstOrder/StiffnessMatrixFirstOrder.h"
#include "../StiffnessMatrix.h"
#include "../../Log/Log.h"
#include "../../Timer/Timer.h"
#include "../../Sparse/Sparse.h"

class StiffnessMatrixSingleCPU : public StiffnessMatrixFirstOrder
{
  //variables
private:
    unsigned concurentThreadsSupported;
    int* simulationPerThread;
    int threadNumber;
  //methods
public:
  StiffnessMatrixSingleCPU(double*, Geometry&, int);
  ~StiffnessMatrixSingleCPU();
  Sparse& GetStiffnessMatrix() override;
  void DisplacementToStrain(double* displacement, double* strain);
private:
};
#endif
