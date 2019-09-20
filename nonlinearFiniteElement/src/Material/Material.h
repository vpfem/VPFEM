#ifndef MATERIAL_H
#define MATERIAL_H

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "Log/Log.h"

class Material {
public:
    Material();
    virtual ~Material();
    virtual void elasticMaterialMatrix(double*, double*) = 0;
    virtual double UpdateMatMatrix(double*, double*, double*) = 0;
};

#endif
