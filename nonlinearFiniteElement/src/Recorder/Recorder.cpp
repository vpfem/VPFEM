#include "Recorder.h"

const void Recorder::matrix(const std::string& fileName , double* ptrArray, unsigned int sizeArray) const {
  std::ofstream outfile;
  std::string fileName2 = std::string("output/") + fileName;
  outfile.open(fileName2);//, iso::out | iso::trunc);
  // I will write a temperary loop it has to change to
  // outfile << youMatrix;
  for (size_t i = 0; i < sizeArray; i++)
    outfile << ptrArray[i] << "\n";
  outfile.close();
}

const void Recorder::SparseMatrix(const std::string& fileName , const Sparse& sp) const {
  std::ofstream outfile;
  std::string fileName2 = std::string("output/") + fileName;
  outfile.open(fileName2);//, iso::out | iso::trunc);
  outfile << sp;
  outfile.close();
}
/*
const void Recorder::SparseMatrix(const std::string& fileName , const Assembly& sp) const {
  std::ofstream outfile;
  std::string fileName2 = std::string("output/") + fileName;
  outfile.open(fileName2);//, iso::out | iso::trunc);
  outfile << sp;
  outfile.close();
}
*/
Recorder& Recorder::File() {
  static Recorder instanceOfRecorder;
  return instanceOfRecorder;
}
