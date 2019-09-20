#include "Application.h"

static double mm = 1.0;//1e-3;
static double MPa = 1.0;//1e6;
//static double kN = 1e3;

void Application::ElasticRun()
{
    //StiffnessMatrixFirstOrder* stiffMat    = new StiffnessMatrixSingleCPU(d_matrix, *cantilever,1);
    StiffnessMatrixFirstOrder* stiffMat    = new StiffnessMatrixGPU(d_matrix, *cantilever,2);
    Sparse &k = stiffMat->GetStiffnessMatrix();
    k.ThrustAssemble2(0.001); // assemble in CPU
    //k.ThrustAssemble3(0.001); // assemble in CPU
    Solver solver(k.get_numberOfRows(), k.get_valueSize(), k.get_value(),
                    (int*) k.get_i(), (int*) k.get_j(),
                    displacement, cantilever->get_Load().get_vector());
    //timer->Start();
    solver.IterativeCgDevice(); // solve in GPU CG
    //solver.IterativeCgPrecondCholDevice(); // solve in GPU Direct
    //timer->Stop();
    //solver.IterativeCgHost(); // solver in CPU
    stiffMat->DisplacementToStrain(displacement, strain);
    delete stiffMat;
    //delete timer;
}

void Application::NonlinearRun()
{
    double error = 1.0;
    double TOL = 1e-5;
    int iterator = 0;
    int MAX_ITERATOR = 100;

    while (error > TOL && iterator < MAX_ITERATOR)
    {
        ElasticRun();
        error = UpdateMaterial();
        //std::cout<<"error in iteration " << iterator << " is: " << error << std::endl;
        iterator++;
    }
}

Application::Application(int e_X, int e_Y)
: elementX(e_X), elementY(e_Y)
{
    Log::Logger().setLevel(Log::LevelInfo);
    //timer = new Timer();
}

Application::~Application()
{
    delete cantilever;
    delete mat;
    cudaFree(displacement);
    cudaFree(d_matrix);
}

void Application::MaterialBuilder()
{
    mat = new Concrete01(21.8*MPa, 1.54*MPa, 100000*MPa, -0.0018, 0.3, 200000*MPa,0,0,0,0,0.02195, 402*MPa,0);
    cudaMallocManaged(&d_matrix, cantilever->get_numberOfElementsG()*6*sizeof(double));
    cudaMallocManaged(&secant_stiffness, cantilever->get_numberOfElementsG()*4*sizeof(double));
    AssignMaterial();
}


void Application::AssignMaterial()
{
    for (int elementNumber = 0; elementNumber < cantilever->get_numberOfElementsG(); elementNumber++)
    {
        double* ptr_element = d_matrix + elementNumber*6; // pointer to the element material matrix 
        double* ptr_secant_stiffness = secant_stiffness + elementNumber*4; // pointer to the secont modules 
        mat->elasticMaterialMatrix(ptr_element, ptr_secant_stiffness);
    }
}

double Application::UpdateMaterial()
{
    double _error = 0.00001;
    for (int elementNumber = 0; elementNumber < cantilever->get_numberOfElementsG(); elementNumber++)
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
    double dimensionX = 4000.0*mm;
    double dimensionY = 400.0*mm;
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

    //int dummyCounter = 0;
    for (int j = 0; j <= ElementNumberY; j++)
    {
        for (int i = 0; i <= ElementNumberX; i++)
        {
            cantilever->node(i*incrementX,j*incrementY);
            //std::cout << "node("<<dummyCounter++ <<") = " << i*incrementX << ", "<< j*incrementY << std::endl;
        }
    }
    cantilever->set_thickness(300.0*mm);

    for (int i = 0; i <= ElementNumberY; i++)
    {
        cantilever->dof->fix(i*(ElementNumberX+1),1,1); // All the nodes in the right side
        //std::cout << "fix(" << i*(ElementNumberX+1) << ")" << std::endl;
    }

    // load 
    double load = 1000;
    double loadVal = load / (ElementNumberY*2);
    int l1 = ((ElementNumberX + 1) * (ElementNumberY + 1))-1; // top right node
    int l2 = ElementNumberX; // bottom right node
    cantilever->load->point(l1,0.0, loadVal);
    cantilever->load->point(l2,0.0, loadVal);
    for (int i = 2; i <= ElementNumberY; i++)
    {
        cantilever->load->point(i*(ElementNumberX+1) - 1 ,0.0, 2*loadVal);
        //std::cout << "load(" << i*(ElementNumberX+1) - 1 << "  , " << 2*loadVal <<" ) " << std::endl;
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
    // dummyCounter = 0;
    for (int i = 0; i < ElementNumberY; i++)
    {
        for (int j = 0; j < ElementNumberX; j++)
        {
            cantilever->meshQuadrilateral(i*xElement+j,i*xElement+1+j,(i+1)*xElement+1+j,(i+1)*xElement+j); 
            //std::cout << "meshQuadrilateral("<<dummyCounter++ <<") = " << i*xElement+j << ", " \
            << i*xElement+1+j << ", "<< (i+1)*xElement+1+j << ", " << (i+1)*xElement+j << std::endl;
        }
    }
    cantilever->modelBuild();


    // allocate memory for displacement and strain
    cudaMallocManaged(&displacement, cantilever->get_Dof().get_freeSize()*sizeof(double));
    int strainPerElement = 3;
    cudaMallocManaged(&strain, strainPerElement*cantilever->get_numberOfElementsG()*sizeof(double));
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

int Application::get_numberOFElements()
{
    return cantilever->get_numberOfElementsG();
}

Geometry* Application::get_Geometry()
{
    return cantilever;
}
