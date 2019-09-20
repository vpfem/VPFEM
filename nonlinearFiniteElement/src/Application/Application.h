#pragma once

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

#include "Timer/Timer.h"
#include "Log/Log.h"

#include "Sparse/Sparse.h"
#include "Geometry/Geometry.h"
#include "Material/Material.h"
#include "Material/ElasticPlainStrain.h"
#include "Material/ElasticPlainStress.h"
#include "Material/Concrete01.h"
#include "StiffnessMatrix/StiffnessMatrix.h"
#include "StiffnessMatrix/StiffnessMatrixFirstOrder/StiffnessMatrixFirstOrder.h"
#include "StiffnessMatrix/StiffnessMatrixParallelCPU/StiffnessMatrixParallelCPU.h"
#include "StiffnessMatrix/StiffnessMatrixSingleCPU/StiffnessMatrixSingleCPU.h"
#include "StiffnessMatrix/StiffnessMatrixGPU/StiffnessMatrixGPU.h"
#include "Recorder/Recorder.h"
#include "Assembly/Assembly.h"
#include "Assembly/AssemblySingleCpu.h"
#include "Assembly/AssemblyParCpu.h"
#include "Solver/Solver.h"


class Application {
private:
    Geometry* cantilever;
    Concrete01* mat;
    double* d_matrix;
    double* displacement;
    double* strain;
    double* secant_stiffness;
    int elementX;
    int elementY;
    Timer* timer;
public:
    Application(int, int);
    ~Application();
    void ElasticRun();
    void ModelBuilder();
    void MaterialBuilder();
    void NonlinearRun();
    double* get_d_matrix();
    double* get_strain();
    double* get_displacement();
    size_t get_numberOFElements();
    Geometry* get_Geometry();
private:
    void AssignMaterial();
    double UpdateMaterial();
};

/*
a[0]= -31150
a[1]= 0
a[2]= 32707.5
a[3]= -96565
a[4]= -31150
a[5]= 0
a[6]= 0
a[7]= 96565
a[8]= 31150
a[9]= -32707.5
a[10]= -15575
a[11]= 31150
a[12]= 0
a[13]= 63857.5
a[14]= 15575
14 % echo END
*/