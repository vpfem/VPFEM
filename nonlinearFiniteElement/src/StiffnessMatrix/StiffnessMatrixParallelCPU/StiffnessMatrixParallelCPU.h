#ifndef STIFFNESSMATRIXPARALLELCPU
#define STIFFNESSMATRIXPARALLELCPU

#include <iostream>
#include <thread>
#include <string>
#include "../StiffnessMatrixFirstOrder/StiffnessMatrixFirstOrder.h"
#include "../StiffnessMatrix.h"
#include "../../Log/Log.h"
#include "../../Timer/Timer.h"
#include "../../Sparse/Sparse.h"

class StiffnessMatrixParallelCPU : public StiffnessMatrixFirstOrder
{
  //variables
private:
  unsigned concurentThreadsSupported;
    int* simulationPerThread;
    int threadNumber;
  //methods
public:
  StiffnessMatrixParallelCPU(double*, Geometry&,   int);
  StiffnessMatrixParallelCPU(double*, Geometry&,   int,   int); // takes number of the desired cores
  ~StiffnessMatrixParallelCPU();
  Sparse& GetStiffnessMatrix() override;
  void DisplacementToStrain(double* displacement, double* strain);
private:  
  void GetStiffnessMatrixForEachThread(  int);
};
#endif
