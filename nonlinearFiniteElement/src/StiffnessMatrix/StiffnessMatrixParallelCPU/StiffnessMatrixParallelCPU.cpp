#include "StiffnessMatrixParallelCPU.h"

StiffnessMatrixParallelCPU::StiffnessMatrixParallelCPU(double* mat, Geometry &geo, int n)
  : StiffnessMatrixParallelCPU(mat, geo, n,(std::thread::hardware_concurrency()-1))
{ 
};

StiffnessMatrixParallelCPU::StiffnessMatrixParallelCPU(double* mat, Geometry &geo, int n, int numberOfCores)
  : StiffnessMatrixFirstOrder(mat, geo, n), concurentThreadsSupported(numberOfCores)
{
  INFO("StiffnessMatrixParallelCPU created with " + std::to_string(concurentThreadsSupported) + " threads");
  simulationPerThread = new int[concurentThreadsSupported];
  for (int i = 0; i < concurentThreadsSupported; i++)
    {
      if (i == concurentThreadsSupported-1)
        simulationPerThread[i] = simulationSize/concurentThreadsSupported + simulationSize%concurentThreadsSupported;
      else
        simulationPerThread[i] = simulationSize/concurentThreadsSupported;
    };
}

StiffnessMatrixParallelCPU::~StiffnessMatrixParallelCPU()
{
  INFO("StiffnessMatrixParallelCPU Deleted");
  delete[] simulationPerThread;

}

Sparse& StiffnessMatrixParallelCPU::GetStiffnessMatrix()
{
  for (int i = 0; i<numberOfElements; i++)
    constantCreator(i, c, geometry->get_x(), geometry->get_y(), geometry->get_mesh());
  //std::cout << "Time spend in CPU using " << concurentThreadsSupported << " cores is: ";
  //Timer timer("");
  std::thread t[concurentThreadsSupported]; // number of threads being used in this program
  for (int i = 0; i<concurentThreadsSupported; i++)
    {
      t[i] = std::thread(&StiffnessMatrixParallelCPU::GetStiffnessMatrixForEachThread,this,i);
    }
  for (int i = 0; i<concurentThreadsSupported; i++)
    {
      t[i].join();
    }
  //stiffMat = assembler(stiffMat, geometry->get_freeDofs(), geometry->get_freeDofs_size());
  return *stiffMat;
}

void StiffnessMatrixParallelCPU::GetStiffnessMatrixForEachThread(int threadId)
{
int counter = threadId*(simulationSize/concurentThreadsSupported);
  for (int i = 0; i<simulationPerThread[threadId]; i++)
    stiffnessMatrixCalculation(counter+i, nipSquared,integrationNode, integrationPos, integrationWeight, c, material, geometry->get_mesh(), stiffMat->value, stiffMat->i, stiffMat->j, geometry->get_Dof().get_free(), geometry->get_thickness());
}

void StiffnessMatrixParallelCPU::DisplacementToStrain(double* displacement, double* strain)
{
    std::thread t[concurentThreadsSupported];
    for (int i = 0; i < numberOfElements; i++)
    {
        DisplacementToStrainEachElement(i, displacement, geometry->get_Dof().get_free(), geometry->get_mesh(), strain, 3, nipSquared,integrationNode, integrationPos, c);
    }
} 