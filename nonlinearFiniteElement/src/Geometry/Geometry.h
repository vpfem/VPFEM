#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cuda.h>
#include <cuda_runtime.h>
#include "../Sparse/Sparse.h"
#include "../Log/Log.h"

// Load Struct
struct Load {
  // variables
  std::vector<double> LoadVector_value;
  std::vector<int> LoadVector_dof_i;
  double* loadVec; 
  // methodes
  Load();
  ~Load();
  void point(int, double, double);
  void line(int, int, double, double);
  void build(int*, int);
  double* get_vector() const;
};

// dof Struct
struct Dof {
  // variables
  int* free; // if zero means it is fixed if numbered is the new dof value
  int freeSize;
  int fixSize;
  std::vector<int> dofFixTemp;
  // methods
  Dof();
  ~Dof();
  void fix(int, int, int);
  void build(int);
  int* get_free() const;
  int get_freeSize() const;
  int get_fixSize () const;
};

// Geometry class
class Geometry {
private:
  int counter;
  std::vector<int> meshTemp;
  std::vector<double> xDim;
  std::vector<double> yDim;
  double thickness;
  double* thickness_array_d;
  int numberOfNodes;
  int numberOfElementsG;
  double* x; double* x_d;
  double* y; double* y_d;
  int* mesh;int* mesh_d;
public:
  Load* load;
  Dof* dof;

public:
  Geometry();
  ~Geometry();
  void node(double, double);
  void modelBuild();
  void meshQuadrilateral(int,int,int,int);
  double* get_x();
  double* get_y();
  int get_x_y_size();
  int get_numberOfElementsG();
  int get_mesh_Size();
  int* get_mesh();
  const Dof& get_Dof() const;
  const Load& get_Load() const;
  const double* get_thickness() const;
  void set_thickness(double t); 
};

#endif
