#include "StiffnessMatrixFirstOrder.h"

void StiffnessMatrixFirstOrder::constantCreator(  int numberElement, double* c, double* x, double* y, int* mesh)
{
    int i = numberElement*6;
    c[i++] = (x[mesh[numberElement*4+0]] - x[mesh[numberElement*4+1]] + x[mesh[numberElement*4+2]] - x[mesh[numberElement*4+3]])/4; 
    c[i++] = (x[mesh[numberElement*4+0]] - x[mesh[numberElement*4+1]] - x[mesh[numberElement*4+2]] + x[mesh[numberElement*4+3]])/4;
    c[i++] = (x[mesh[numberElement*4+0]] - x[mesh[numberElement*4+3]] - x[mesh[numberElement*4+2]] + x[mesh[numberElement*4+1]])/4;
    c[i++] = (y[mesh[numberElement*4+0]] - y[mesh[numberElement*4+1]] + y[mesh[numberElement*4+2]] - y[mesh[numberElement*4+3]])/4;
    c[i++] = (y[mesh[numberElement*4+0]] - y[mesh[numberElement*4+1]] - y[mesh[numberElement*4+2]] + y[mesh[numberElement*4+3]])/4;
    c[i++] = (y[mesh[numberElement*4+0]] - y[mesh[numberElement*4+3]] - y[mesh[numberElement*4+2]] + y[mesh[numberElement*4+1]])/4;
    // defined the constants c1x to c3y
};

void StiffnessMatrixFirstOrder::stiffnessMatrixCalculation(  int numberElement, 
    int nip ,double* in, int* ip, double* iw, double* c, double* D, 
    int* mesh, double* k, int* i_index, int* j_index, 
    int* dofFree, const double* thickness)
// numberElement -> the element number needed to be calculated
// nip is the number of integration point squared.
// in is the integrationNode
// ip -> integrationPos
// iw -> integrationWeight
// c -> constants
// mesh -> a pointer to mesh array 
// D -> material matrix
// k -> pointer to stiffness matrix
// i_index -> pointer to DOFi
// j_index -> pointer to DOFj
// dofFree -> lists the free dofs and value for new Dofs
// thickness -> thickness of each element
{
    // define local stiffness Matrix
    double kLocal[64] = {}; 
    for (  int noIP = 0; noIP < nip; noIP++) { 
        // noIP -> the integration point number needed to be calculated
        double XI = in[ip[2*noIP]]; double YI = in[ip[2*noIP+1]];
        // Jacobian
        double J11 = c[numberElement*6+0]*YI-c[numberElement*6+1]; double J12 = c[numberElement*6+3]*YI-c[numberElement*6+4];
        double J21 = c[numberElement*6+0]*XI-c[numberElement*6+2]; double J22 = c[numberElement*6+3]*XI-c[numberElement*6+5];
        double detJ = J11*J22-J12*J21;
        // thickness of the element
        double t = thickness[numberElement];
        double WeightPerDetJ = t*(iw[ip[2*noIP]]*iw[ip[2*noIP+1]])/detJ;
        // derveativs of the shape function N1x N2x ... N1y N2y ...
        double Ni[8] = {J22*( YI-1)/4 -  J12*( XI-1)/4, J22*(-YI+1)/4 -  J12*(-XI-1)/4, \
                        J22*( YI+1)/4 -  J12*( XI+1)/4, J22*(-YI-1)/4 -  J12*(-XI+1)/4, \
                        J11*( XI-1)/4 -  J21*( YI-1)/4, J11*(-XI-1)/4 -  J21*(-YI+1)/4, \
                        J11*( XI+1)/4 -  J21*( YI+1)/4, J11*(-XI+1)/4 -  J21*(-YI-1)/4};
        // multiplication of shape functions N1x^2 N1x*N2x ....
        double N[64];
        int counterN = 0;
        for (  int i = 0; i < 8; i++) 
        {
            for (  int j = i; j < 8 ; j++)
                N[counterN++] = Ni[i]*Ni[j];
        };
        // find the position to start filling the stiffness matrix
        // writes all 64 components of the 8 by 8 symmteric stiffness Matrix (storing upper triangle)
        int entry = 0;
                 kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[0] + 2*D[numberElement*6+4]*N[4] + D[numberElement*6+2]*N[26]); //1,1             
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[0] + D[numberElement*6+5]*N[26] + D[numberElement*6+3]*N[4] + D[numberElement*6+2]*N[4]); //2,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[1] + D[numberElement*6+4]*N[5] + D[numberElement*6+4]*N[11] + D[numberElement*6+2]*N[27]);//3,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[1] + D[numberElement*6+3]*N[5] + D[numberElement*6+2]*N[11] + D[numberElement*6+5]*N[27]);//4,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[2] + D[numberElement*6+4]*N[6] + D[numberElement*6+4]*N[17] + D[numberElement*6+2]*N[28]);//5,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[2] + D[numberElement*6+3]*N[6] + D[numberElement*6+2]*N[17] + D[numberElement*6+5]*N[28]);//6,1
        entry++; kLocal[entry] = kLocal[entry]  + WeightPerDetJ*(D[numberElement*6+0]*N[3] + D[numberElement*6+4]*N[7] + D[numberElement*6+4]*N[22] + D[numberElement*6+2]*N[29]); //7,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[3] + D[numberElement*6+3]*N[7] + D[numberElement*6+2]*N[22] + D[numberElement*6+5]*N[29]);//8,1

        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[0] + D[numberElement*6+5]*N[26] + D[numberElement*6+3]*N[4] + D[numberElement*6+2]*N[4]); //2,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[0] + 2*D[numberElement*6+5]*N[4] + D[numberElement*6+1]*N[26]);           //2,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[1] + D[numberElement*6+3]*N[11] + D[numberElement*6+2]*N[5] + D[numberElement*6+5]*N[27]);//3,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[1] + D[numberElement*6+5]*N[5] + D[numberElement*6+5]*N[11] + D[numberElement*6+1]*N[27]);//4,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[2] + D[numberElement*6+3]*N[17] + D[numberElement*6+2]*N[6] + D[numberElement*6+5]*N[28]);//5,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[2] + D[numberElement*6+5]*N[6] + D[numberElement*6+5]*N[17] + D[numberElement*6+1]*N[28]);//6,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[3] + D[numberElement*6+3]*N[22] + D[numberElement*6+2]*N[7] + D[numberElement*6+5]*N[29]); //7,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[3] + D[numberElement*6+5]*N[7] + D[numberElement*6+5]*N[22] + D[numberElement*6+1]*N[29]);//8,2

        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[1] + D[numberElement*6+4]*N[5] + D[numberElement*6+4]*N[11] + D[numberElement*6+2]*N[27]);//3,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[1] + D[numberElement*6+3]*N[11] + D[numberElement*6+2]*N[5] + D[numberElement*6+5]*N[27]);//3,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[8] + 2*D[numberElement*6+4]*N[12] + D[numberElement*6+2]*N[30]);          //3,3
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[8] + D[numberElement*6+5]*N[30] + D[numberElement*6+3]*N[12] + D[numberElement*6+2]*N[12]);//4,3
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[9] + D[numberElement*6+4]*N[13] + D[numberElement*6+4]*N[18] + D[numberElement*6+2]*N[31]);//5,3
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[9] + D[numberElement*6+3]*N[13] + D[numberElement*6+2]*N[18] + D[numberElement*6+5]*N[31]);//6,3
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[10] + D[numberElement*6+4]*N[14] + D[numberElement*6+4]*N[23] + D[numberElement*6+2]*N[32]);//7,3
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[10] + D[numberElement*6+3]*N[14] + D[numberElement*6+2]*N[23] + D[numberElement*6+5]*N[32]);//8,3

        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[1] + D[numberElement*6+3]*N[5] + D[numberElement*6+2]*N[11] + D[numberElement*6+5]*N[27]);//4,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[1] + D[numberElement*6+5]*N[5] + D[numberElement*6+5]*N[11] + D[numberElement*6+1]*N[27]);//4,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[8] + D[numberElement*6+5]*N[30] + D[numberElement*6+3]*N[12] + D[numberElement*6+2]*N[12]);//4,3
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[8] + 2*D[numberElement*6+5]*N[12] + D[numberElement*6+1]*N[30]);           //4,4
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[9] + D[numberElement*6+3]*N[18] + D[numberElement*6+2]*N[13] + D[numberElement*6+5]*N[31]);//5,4
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[9] + D[numberElement*6+5]*N[13] + D[numberElement*6+5]*N[18] + D[numberElement*6+1]*N[31]);//6,4
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[10] + D[numberElement*6+3]*N[23] + D[numberElement*6+2]*N[14] + D[numberElement*6+5]*N[32]);//7,4
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[10] + D[numberElement*6+5]*N[14] + D[numberElement*6+5]*N[23] + D[numberElement*6+1]*N[32]);//8,4

        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[2] + D[numberElement*6+4]*N[6] + D[numberElement*6+4]*N[17] + D[numberElement*6+2]*N[28]);//5,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[2] + D[numberElement*6+3]*N[17] + D[numberElement*6+2]*N[6] + D[numberElement*6+5]*N[28]);//5,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[9] + D[numberElement*6+4]*N[13] + D[numberElement*6+4]*N[18] + D[numberElement*6+2]*N[31]);//5,3
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[9] + D[numberElement*6+3]*N[18] + D[numberElement*6+2]*N[13] + D[numberElement*6+5]*N[31]);//5,4
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[15] + 2*D[numberElement*6+4]*N[19] + D[numberElement*6+2]*N[33]);//5,5
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[15] + D[numberElement*6+5]*N[33] + D[numberElement*6+3]*N[19] + D[numberElement*6+2]*N[19]);//6,5
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[16] + D[numberElement*6+4]*N[20] + D[numberElement*6+4]*N[24] + D[numberElement*6+2]*N[34]);//7,5
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[16] + D[numberElement*6+3]*N[20] + D[numberElement*6+2]*N[24] + D[numberElement*6+5]*N[34]);//8,5

        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[2] + D[numberElement*6+3]*N[6] + D[numberElement*6+2]*N[17] + D[numberElement*6+5]*N[28]);//6,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[2] + D[numberElement*6+5]*N[6] + D[numberElement*6+5]*N[17] + D[numberElement*6+1]*N[28]);//6,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[9] + D[numberElement*6+3]*N[13] + D[numberElement*6+2]*N[18] + D[numberElement*6+5]*N[31]);//6,3
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[9] + D[numberElement*6+5]*N[13] + D[numberElement*6+5]*N[18] + D[numberElement*6+1]*N[31]);//6,4
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[15] + D[numberElement*6+5]*N[33] + D[numberElement*6+3]*N[19] + D[numberElement*6+2]*N[19]);//6,5
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[15] + 2*D[numberElement*6+5]*N[19] + D[numberElement*6+1]*N[33]);//6,6
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[16] + D[numberElement*6+3]*N[24] + D[numberElement*6+2]*N[20] + D[numberElement*6+5]*N[34]);//7,6
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[16] + D[numberElement*6+5]*N[20] + D[numberElement*6+5]*N[24] + D[numberElement*6+1]*N[34]);//8,6

        entry++; kLocal[entry] = kLocal[entry]  + WeightPerDetJ*(D[numberElement*6+0]*N[3] + D[numberElement*6+4]*N[7] + D[numberElement*6+4]*N[22] + D[numberElement*6+2]*N[29]); //7,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[3] + D[numberElement*6+3]*N[22] + D[numberElement*6+2]*N[7] + D[numberElement*6+5]*N[29]); //7,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[10] + D[numberElement*6+4]*N[14] + D[numberElement*6+4]*N[23] + D[numberElement*6+2]*N[32]);//7,3
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[10] + D[numberElement*6+3]*N[23] + D[numberElement*6+2]*N[14] + D[numberElement*6+5]*N[32]);//7,4
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[16] + D[numberElement*6+4]*N[20] + D[numberElement*6+4]*N[24] + D[numberElement*6+2]*N[34]);//7,5
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[16] + D[numberElement*6+3]*N[24] + D[numberElement*6+2]*N[20] + D[numberElement*6+5]*N[34]);//7,6
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+0]*N[21] + 2*D[numberElement*6+4]*N[25] + D[numberElement*6+2]*N[35]);//7,7
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[21] + D[numberElement*6+5]*N[35] + D[numberElement*6+3]*N[25] + D[numberElement*6+2]*N[25]);//8,7

        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[3] + D[numberElement*6+3]*N[7] + D[numberElement*6+2]*N[22] + D[numberElement*6+5]*N[29]);//8,1
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[3] + D[numberElement*6+5]*N[7] + D[numberElement*6+5]*N[22] + D[numberElement*6+1]*N[29]);//8,2
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[10] + D[numberElement*6+3]*N[14] + D[numberElement*6+2]*N[23] + D[numberElement*6+5]*N[32]);//8,3
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[10] + D[numberElement*6+5]*N[14] + D[numberElement*6+5]*N[23] + D[numberElement*6+1]*N[32]);//8,4
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[16] + D[numberElement*6+3]*N[20] + D[numberElement*6+2]*N[24] + D[numberElement*6+5]*N[34]);//8,5
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[16] + D[numberElement*6+5]*N[20] + D[numberElement*6+5]*N[24] + D[numberElement*6+1]*N[34]);//8,6
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+4]*N[21] + D[numberElement*6+5]*N[35] + D[numberElement*6+3]*N[25] + D[numberElement*6+2]*N[25]);//8,7
        entry++; kLocal[entry] = kLocal[entry] + WeightPerDetJ*(D[numberElement*6+2]*N[21] + 2*D[numberElement*6+5]*N[25] + D[numberElement*6+1]*N[35]);//8,8
    }
    int counter = 64*(numberElement);
    int kcounter = 0, dof_i,dof_j;
    for (  int i = 0; i < 8; i++) {
        for (  int j = 0; j < 8; j++) 
        {
            dof_i = (mesh[numberElement*4+i/2])*2+i%2+1; // this (mesh[numberElement*4+i/2])*2+i%2+1 the required dof however we need to look at dofFree which the index starts at zero
            dof_j = (mesh[numberElement*4+j/2])*2+j%2+1;
            if (dofFree[dof_i-1] && dofFree[dof_j-1]) 
            {
                i_index[counter] = dofFree[dof_i-1];
                j_index[counter] = dofFree[dof_j-1];
                k[counter] = kLocal[kcounter];
            }
            //rowPtr[i_index[counter]] = rowPtr[i_index[counter]] + 1;
            counter++;
            kcounter++;
        }
    }
}

StiffnessMatrixFirstOrder::StiffnessMatrixFirstOrder(double* mat, Geometry& geo, int n)
    :StiffnessMatrix(mat,geo,n)
{
    INFO("StiffnessMatrixFirstOrder Created by CPU");
    sizeStiffMatPerEle = 64;
    stiffMatSize = numberOfElements*sizeStiffMatPerEle;
    simulationSize = numberOfElements;
    stiffMat = new Sparse(stiffMatSize,geometry->get_Dof().get_freeSize());
    cudaMallocManaged(&c,numberOfElements*6*sizeof(double));
};

StiffnessMatrixFirstOrder::~StiffnessMatrixFirstOrder() {
    INFO("StiffnessMatrixFirstOrder Deleted by CPU");
    cudaFree(c);
    //delete stiffMat;
}

int StiffnessMatrixFirstOrder::GetStiffnessMatrixSize()
{
    return stiffMatSize;
}


void StiffnessMatrixFirstOrder::DisplacementToStrainEachElement(int elementNumber,\
    double * displacement, int * dofFree, int * mesh, double* strain,\
  int strainPerElement, int nip ,double* in, int* ip, double* c)
// elementNumber
// displacement: pointer to the displacement of all free nodes
// dofFree: int if fix entry is 0
// mesh: mesh diagram which is an int *
// nip is the number of integration point squared.
// in: integration node
// ip: integration position 
// c: constants calculated using x and y
{
    // make an array of 8 holding the displacements of all nodes in the element 
    double Node;
    double d[8] = {0};
  int counter = 0;
    for (int i = 0; i < 4; i++)
    {
        Node = mesh[elementNumber*4 + i];
        for (int j = 2*Node+1; j < 2*Node+3; j++)
        {
            if (dofFree[j-1])
            {
                d[counter] = displacement[dofFree[j-1]-1];
            }
            counter++;
        }
    }
    //ARRAY(d,8);
    // the best way to do this is epsilon = B*d
    double Ni[8] = {0};
    for (  int noIP = 0; noIP < nip; noIP++) { 
        // noIP -> the integration point number needed to be calculated
        double XI = in[ip[2*noIP]]; double YI = in[ip[2*noIP+1]];
        // Jacobian
        double J11 = c[elementNumber*6+0]*YI-c[elementNumber*6+1]; double J12 = c[elementNumber*6+3]*YI-c[elementNumber*6+4];
        double J21 = c[elementNumber*6+0]*XI-c[elementNumber*6+2]; double J22 = c[elementNumber*6+3]*XI-c[elementNumber*6+5];
        double detJ = J11*J22-J12*J21;

        // derveativs of the shape function N1x N2x ... N1y N2y ...
        Ni[0] = Ni[0] + (J22*( YI-1)/4 -  J12*( XI-1)/4)/(detJ*nip);
        Ni[1] = Ni[1] + (J22*(-YI+1)/4 -  J12*(-XI-1)/4)/(detJ*nip);
        Ni[2] = Ni[2] + (J22*( YI+1)/4 -  J12*( XI+1)/4)/(detJ*nip);
        Ni[3] = Ni[3] + (J22*(-YI-1)/4 -  J12*(-XI+1)/4)/(detJ*nip);
        Ni[4] = Ni[4] + (J11*( XI-1)/4 -  J21*( YI-1)/4)/(detJ*nip);
        Ni[5] = Ni[5] + (J11*(-XI-1)/4 -  J21*(-YI+1)/4)/(detJ*nip);
        Ni[6] = Ni[6] + (J11*( XI+1)/4 -  J21*( YI+1)/4)/(detJ*nip);
        Ni[7] = Ni[7] + (J11*(-XI+1)/4 -  J21*(-YI-1)/4)/(detJ*nip);
    }
    // this is where we perform B*d
    strain[elementNumber*3 + 0] = d[0]*Ni[0] + d[2]*Ni[1]  + d[4]*Ni[2] + d[6]*Ni[3];
    strain[elementNumber*3 + 1] = d[1]*Ni[4] + d[3]*Ni[5]  + d[5]*Ni[6] + d[7]*Ni[7];
    strain[elementNumber*3 + 2] = d[0]*Ni[4] + d[2]*Ni[5]  + d[4]*Ni[6] + d[6]*Ni[7]+ \
                                  d[1]*Ni[0] + d[3]*Ni[1]  + d[5]*Ni[2] + d[7]*Ni[3];
    //ARRAY(strain,3);
}







