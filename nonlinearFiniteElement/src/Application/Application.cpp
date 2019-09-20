#include "Application.h"

static double mm = 1.0;//1e-3;
static double MPa = 1.0;//1e6;
//static double kN = 1e3;

void Application::ElasticRun()
{
    //StiffnessMatrixFirstOrder* stiffMat    = new StiffnessMatrixSingleCPU(d_matrix, *cantilever,4);
    StiffnessMatrixFirstOrder* stiffMat    = new StiffnessMatrixGPU(d_matrix, *cantilever,4);
    Sparse &k = stiffMat->GetStiffnessMatrix();
    //k.STLAssemble(0.001); // assemble in CPU
    k.ThrustAssemble2(0.001); // assembly in GPU
    Solver solver(k.get_numberOfRows(), k.get_valueSize(), k.get_value(),
                    (int*) k.get_i(), (int*) k.get_j(),
                    displacement, cantilever->get_Load().get_vector());
    solver.IterativeCgDevice(); // solve in GPU
    //solver.IterativeCgHost(); // solver in CPU
    stiffMat->DisplacementToStrain(displacement, strain);
    delete stiffMat;
}

void Application::NonlinearRun()
{
    // allocate memory for displacement and strain
    cudaMallocManaged(&displacement, cantilever->get_Dof().get_freeSize()*sizeof(double));
    int strainPerElement = 3;
    cudaMallocManaged(&strain, strainPerElement*cantilever->get_numberOfElementsG()*sizeof(double));

    double error = 1.0;
    double TOL = 1e-5;
    int iterator = 0;
    int MAX_ITERATOR = 100;

    while (error > TOL && iterator < MAX_ITERATOR)
    {
        ElasticRun();
        error = UpdateMaterial();
        // printing strain in each step
        //std::cout<<"error in iteration " << iterator << " is: " << error << std::endl;
        iterator++;
    }
}

Application::Application(int e_X, int e_Y)
: elementX(e_X), elementY(e_Y)
{
    Log::Logger().setLevel(Log::LevelInfo);
    //timer = new Timer("App timer: ");
}

Application::~Application()
{
    delete cantilever;
    delete mat;
    //delete timer;
    cudaFree(displacement);
    cudaFree(d_matrix);
    //cudaFree(displacement);
    //cudaFree(strain);
}

void Application::MaterialBuilder()
{
    mat = new Concrete01(21.8*MPa, 1.54*MPa, 24200*MPa, -0.0018, 0.3, 200000*MPa,0,0,0,0,0.02195, 402*MPa,0);
    cudaMallocManaged(&d_matrix, cantilever->get_numberOfElementsG()*6*sizeof(double));
    cudaMallocManaged(&secant_stiffness, cantilever->get_numberOfElementsG()*4*sizeof(double));
    AssignMaterial();
}


void Application::AssignMaterial()
{
    for (size_t elementNumber = 0; elementNumber < cantilever->get_numberOfElementsG(); elementNumber++)
    {
        double* ptr_element = d_matrix + elementNumber*6; // pointer to the element material matrix 
        double* ptr_secant_stiffness = secant_stiffness + elementNumber*4; // pointer to the secont modules 
        mat->elasticMaterialMatrix(ptr_element, ptr_secant_stiffness);
    }
}

double Application::UpdateMaterial()
{
    double _error = 0.00001;
    for (size_t elementNumber = 0; elementNumber < cantilever->get_numberOfElementsG(); elementNumber++)
    {
        double* ptr_element = d_matrix + elementNumber*6; // pointer to the element material matrix 
        double* ptr_secant_stiffness = secant_stiffness + elementNumber*4; // pointer to the secont modules
        double* ptr_strain = strain + elementNumber*3; // pointer to the element strain vector 
        double error = mat->UpdateMatMatrix(ptr_element, ptr_secant_stiffness, ptr_strain);
        //std::cout<<"error: " << error << std::endl;
        _error = (error > _error) ? error : _error;
    }
    return _error;
}

void Application::ModelBuilder()
{
    // define units
    // GEOMETRY (building the cantilever)
    double dimensionX = 890.0*mm;
    double dimensionY = 890.0*mm;
    cantilever = new Geometry();
        /*
        2 --- 3 
        |  1  | 
        0 --- 1 
    */
    int ElementNumberX = elementX;
    int ElementNumberY = elementY;
    double incrementX = dimensionX/ElementNumberX;
    double incrementY = dimensionY/ElementNumberY;

    for (int j = 0; j <= ElementNumberY; j++)
    {
        for (int i = 0; i <= ElementNumberX; i++)
        {
            cantilever->node(i*incrementX,j*incrementY);
            //std::cout << "node("<<dummyCounter++ <<") = " << i*incrementX << ", "<< j*incrementY << std::endl;
        }
    }
    cantilever->set_thickness(70.0*mm);

    cantilever->dof->fix(0,1,1); // bottom right node
    cantilever->dof->fix(ElementNumberX,0,1); // bottom left node
    
    double p1 = 1.0*70.0; // Distributed shear load
    double p2 = 3.1*70.0; // Distributed normal load
    
    int nodeX = ElementNumberX+1;
    int nodeY = ElementNumberY+1;
    int incr = 1;
    int n = 1;
    double xLoad, yLoad;
    xLoad = - p1*incrementX/2;
    yLoad = 0.0;

    for (int i = 1; i < 2*(nodeX+nodeY)-3; i++)
    {
        cantilever->load->point(n-incr,xLoad, yLoad);
        cantilever->load->point(n,xLoad, yLoad);
        if (i == nodeX-1)
        {
            incr = nodeX;
            xLoad = + p2*incrementY/2;
            yLoad = + p1*incrementY/2;
        } else if (i == nodeX-1 + nodeY-1)
        {
            incr = -1;
            xLoad = + p1*incrementX/2;
            yLoad = 0.0;
        } else if (i == 2*(nodeX-1) + nodeY-1)
        {
            incr = -nodeX;
            xLoad = - p2*incrementY/2;
            yLoad = - p1*incrementY/2;
        }
        n = n + incr;
    }

    // Mesh the Cantilever beam
        /*
        6 --- 7 --- 8
        |  3  |  4  |
        3 --- 4 --- 5
        |  1  |  2  |
        0 --- 1 --- 2
    */
     int xElement = ElementNumberX+1;
    for (int i = 0; i < ElementNumberY; i++)
    {
        for (int j = 0; j < ElementNumberX; j++)
        {
            cantilever->meshQuadrilateral(i*xElement+j,i*xElement+1+j,(i+1)*xElement+1+j,(i+1)*xElement+j); 
            //std::cout << "meshQuadrilateral("<<dummyCounter++ <<") = " << i*xElement+j << ", "
            //<< i*xElement+1+j << ", "<< (i+1)*xElement+1+j << ", " << (i+1)*xElement+j << std::endl;
        }
    }
    cantilever->modelBuild();
}

double* Application::get_d_matrix()
{
    return d_matrix;
}
double* Application::get_strain()
{
    return strain;
}
double* Application::get_displacement()
{
    return displacement;
}

size_t Application::get_numberOFElements()
{
    return cantilever->get_numberOfElementsG();
}

Geometry* Application::get_Geometry()
{
    return cantilever;
}
