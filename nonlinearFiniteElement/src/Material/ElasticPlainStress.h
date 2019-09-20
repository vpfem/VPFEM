#pragma once 

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>
#include "../Log/Log.h"
#include "Material.h"

class ElasticPlainStress: public Material {
private:
  double ElasticModulus; double poissonRatio;

public:
  ElasticPlainStress(double, double);
  ~ElasticPlainStress();
  virtual void elasticMaterialMatrix(double*, double *) override;
  virtual double UpdateMatMatrix(double*, double*, double*) override;
};