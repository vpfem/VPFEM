#ifndef ELEMENT_H
#define ELEMENT_H

#include <iostream>
#include <vector>
#include <thread>
#include "Timer.h"
#include "Material.h"
#include "Geometry.h"
#include "Log.h"
#include "cs.h"


class Element
{
  // variables
private:
  Log* log;
  unsigned int numberOfIntegrationPoint;
  unsigned int n2IP; //number of IP squared
  unsigned int numberOfElements; unsigned int *NOE; //number of elements
  std::vector<double> integrationNode;
  std::vector<double> integrationWeight;
  Material* material;
  Geometry* geometry;
  std::vector<unsigned int> DOF_1; // DOF_1 a vector rows: 36*NOE  
  std::vector<unsigned int> DOF_2; // (DOF_1, DOF_2) this shows the dofs related to stiffness matrix for the project
  std::vector<std::vector<double>> c; // c1x to x3y -> builds vectors of constants required in calculation of stiffness matrix each row 6 constant for each element
  unsigned int listIndex[36];//list of 36 members of stiffness matrix 
public:
  std::vector<std::vector<double>> stiffnessMatrix; // stiffness matrix a 2D vector rows: 36*NOE and columns: NIP^2
  // Methods
private:

public:
  Element(Log&, Material&, Geometry&, unsigned int);
  void integrationPoint(); // creates integration points and weight
  void stiffnessMatrixCreator(unsigned int);
  void DOF_Creator();
  void stiffnessMatrixFirtOrder(std::vector<std::vector<double>>&, unsigned int, unsigned int);
  void constantCreator();
  void listInitializer();
  // Static Methods
private:
public:

  
};
#endif
