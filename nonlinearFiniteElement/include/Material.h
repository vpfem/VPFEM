#ifndef MATERIAL_H
#define MATERIAL_H

#include <iostream>
#include <vector>
#include <string>
#include "Log.h"

// Geometry class                                                               
class Material {
private:
  float ElasticModulus;
  float poissonRatio;
  Log* log;
public:
  std::vector<float> materialMatrix;
  Material(Log&, float, float);
  void elasticMaterial();
};

#endif
