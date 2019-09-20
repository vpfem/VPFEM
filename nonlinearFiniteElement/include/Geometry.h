#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <vector>
#include <string>
#include "Log.h"

// Geometry class
class Geometry {
private:
  float dimentionX_ ;
  float dimentionY_ ;
  unsigned int numberOfElementX_;
  unsigned int numberOfElementY_;
  unsigned int numberOfElementG_;
  unsigned int numberOfNodes_;
  Log* log;
public:
  std::vector<std::vector<unsigned int>> mesh;
  std::vector<float> x;
  std::vector<float> y;
  Geometry(Log&, float, float, unsigned int, unsigned int);
  Geometry(const Geometry& geometry);
  void nodeDimentionVector();
  void meshCantilever();
};

#endif
