#ifndef RECORDER_H
#define RECORDER_H

#include <iostream>
#include <string>
#include <fstream>
#include "../Sparse/Sparse.h"
#include "../Assembly/Assembly.h"

class Recorder {


public:
  const void matrix(const std::string& , double* , unsigned int) const;
  const void SparseMatrix(const std::string& fileName , const Sparse& sp) const;
  const void SparseMatrix(const std::string& fileName , const Assembly& sp) const;
  static Recorder& File();
  
};

#endif 
