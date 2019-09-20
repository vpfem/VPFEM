#include "testBench.h"

static double mm = 1e-3;
static double MPa = 1e6;
static double kN = 1e3;

int main() {
    // -- start logging
    Log::Logger().setLevel(Log::LevelInfo);
    // -- GEOMETRY (building the cantilever)
    Geometry p_cantilever = geometry();
    // -- Build material vector (E, nu)
    //ElasticPlainStress mat(24200*MPa, 0.3);
    // building material vector try to fix it
    Concrete01 mat1(21.8*MPa, 1.54*MPa, 24200*MPa, -0.0018, 0.3, 200000*MPa,0,0,0,0,0.02195, 402*MPa,0);
    double* d_matrix;
    cudaMallocManaged(&d_matrix, p_cantilever.get_numberOfElementsG()*6*sizeof(double));
    for (int i = 0; i < p_cantilever.get_numberOfElementsG(); i++)
        AssignMaterial(i,mat1, d_matrix);
    for (int i = 0; i < 2; i++)
    {
    INFO(d_matrix,6);
    StiffnessMatrixFirstOrder* stiffMat    = new StiffnessMatrixSingleCPU(d_matrix,p_cantilever,2);
    //StiffnessMatrixFirstOrder* stiffMat    = new StiffnessMatrixGPU(mat.materialMatrix,p_cantilever,4);
    // -- Calculate the stiffness matrix and do the assembly
    Sparse &k = stiffMat->GetStiffnessMatrix();
    // -- Assembly
    //Recorder::File().SparseMatrix("k.out", k);
    k.STLAssemble(0.001);
    //k.ThrustAssemble2(0.001);
    //Recorder::File().SparseMatrix("kAss.out", k);
    // -- Solver
    double* displacement;
    cudaMallocManaged(&displacement, k.get_numberOfRows()*sizeof(double));
    Solver(k.get_numberOfRows(), k.get_valueSize(), k.get_value(), (int*) k.get_i(), (int*) k.get_j(),
        displacement, p_cantilever.get_Load().get_vector());
    stiffMat->set_displacement(displacement);
    // -- Recorder
    INFO(displacement[k.get_numberOfRows()-1]);
    Recorder::File().matrix("displacement.out", displacement, k.get_numberOfRows());

    // calculate the strain of the model
    stiffMat->DisplacementToStrain();
    Recorder::File().matrix("strain.out", stiffMat->get_strain(), 3*p_cantilever.get_numberOfElementsG());
    INFO(stiffMat->get_strain(), 3*p_cantilever.get_numberOfElementsG());

    // update the material stiffness matrix
    for (int i = 0; i < p_cantilever.get_numberOfElementsG(); i++)
        updateMaterial(i, mat1, d_matrix, stiffMat->get_strain());

    
    //cudaFree(d_matrix);
    //cudaFree(displacement);
    delete stiffMat;
    }
    return 0;
}

void AssignMaterial(int elementNumber, Material& mat, double* matMatrix)
{
    double* ptr_element = matMatrix + elementNumber*6; // pointer to the element material matrix 
    mat.elasticMaterialMatrix(ptr_element);
}

void updateMaterial(int elementNumber, Material& mat, double* matMatrix, double* strain)
{
    double* ptr_element = matMatrix + elementNumber*6; // pointer to the element material matrix 
    double* ptr_strain = strain + elementNumber*3; // pointer to the element strain vector 
    mat.UpdateMatMatrix(ptr_element, ptr_strain);
}

Geometry& geometry() {
    // define units
    // GEOMETRY (building the cantilever)
    double dimensionX = 890.0*mm;
    double dimensionY = 890.0*mm;
    Geometry* cantilever = new Geometry();
        /*
        2 --- 3 
        |  1  | 
        0 --- 1 
    */
    int ElementNumberX = 1;
    int ElementNumberY = 1;
    double incrementX = dimensionX/ElementNumberX;
    double incrementY = dimensionX/ElementNumberY;
    int dummyCounter = 0;
    for (int j = 0; j <= ElementNumberY; j++)
    {
        for (int i = 0; i <= ElementNumberX; i++)
        {
            cantilever->node(i*incrementX,j*incrementY);
            //std::cout << "node("<<dummyCounter++ <<") = " << i*incrementX << ", "<< j*incrementY << std::endl;
        }
    }

    cantilever->set_thickness(70.0*mm);

    //for (unsigned int i = 0; i < numberOfElementX; i++)
    //  cantilever->dof->fix(i*(numberOfElementY+1),1,1);
    cantilever->dof->fix(0,1,1); // bottom right node
    cantilever->dof->fix(ElementNumberX,0,1); // bottom left node
    
    cantilever->load->point(ElementNumberX,65.41*kN,0.0); // bottom left node
    
    int topLeftNode = ElementNumberX*ElementNumberY + ElementNumberY;
    cantilever->load->point(topLeftNode,-65.41*kN,-31.15*kN); // top left Node
    int topRightNode = ElementNumberX*ElementNumberY + ElementNumberX + ElementNumberY;
    cantilever->load->point(topRightNode,127.71*kN,31.15*kN); // top right Node
    
    // Mesh the Cantilever beam
        /*
        6 --- 7 --- 8
        |  3  |  4  |
        3 --- 4 --- 5
        |  1  |  2  |
        0 --- 1 --- 2
    */
     int xElement = ElementNumberX+1;
    for (unsigned int i = 0; i < ElementNumberY; i++)
    {
        for (unsigned int j = 0; j < ElementNumberX; j++)
        {
            cantilever->meshQuadrilateral(i*xElement+j,i*xElement+1+j,(i+1)*xElement+1+j,(i+1)*xElement+j); 
            //std::cout << "meshQuadrilateral("<<dummyCounter++ <<") = " << i*xElement+j << ", "
            //<< i*xElement+1+j << ", "<< (i+1)*xElement+1+j << ", " << (i+1)*xElement+j << std::endl;
        }
    }
    cantilever->modelBuild();
    return *cantilever;
}

