#include "StiffnessMatrixSingleCPU.h"

StiffnessMatrixSingleCPU::StiffnessMatrixSingleCPU(double* mat, Geometry &geo, int n)
    : StiffnessMatrixFirstOrder(mat, geo, n)
{ 
    INFO("StiffnessMatrixSingleCPU created");
};

StiffnessMatrixSingleCPU::~StiffnessMatrixSingleCPU()
{
     INFO("StiffnessMatrixSingleCPU Deleted");
}

Sparse& StiffnessMatrixSingleCPU::GetStiffnessMatrix()
{
    for (int i = 0; i<numberOfElements; i++)
        constantCreator(i, c, geometry->get_x(), geometry->get_y(), geometry->get_mesh());
    //Timer timer("Time spend in CPU using single core: ");
    for (int i = 0; i<numberOfElements; i++)
        stiffnessMatrixCalculation(i, nipSquared,integrationNode, integrationPos, integrationWeight, c, material, geometry->get_mesh(), stiffMat->value, stiffMat->i, stiffMat->j, geometry->get_Dof().get_free(),geometry->get_thickness());
    return *stiffMat;
}

void StiffnessMatrixSingleCPU::DisplacementToStrain(double* displacement, double* strain)
{
    for (int i = 0; i < numberOfElements; i++)
    {
        DisplacementToStrainEachElement(i, displacement, geometry->get_Dof().get_free(), geometry->get_mesh(), strain, 3, nipSquared,integrationNode, integrationPos, c);
    }
} 