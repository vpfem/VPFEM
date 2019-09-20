#pragma once 

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "../Log/Log.h"
#include "Material.h"

class ElasticPlainStrain: public Material {
private:
    double ElasticModulus; double poissonRatio;
public:
    ElasticPlainStrain(double, double);
    ~ElasticPlainStrain();
    virtual void elasticMaterialMatrix(double*, double*) override;
};