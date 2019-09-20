#ifndef ASSEMBLYPARCPU_H
#define ASSEMBLYPARCPU_H

#include "AssemblySingleCpu.h"
#include <algorithm>
#include <thread>
#include <cmath>

class Pair 
{
public:
    int x;
    int y;
    Pair();
    Pair(int, int);
};
std::ostream& operator<< (std::ostream &, Pair const& );


class AssemblyParCpu : public AssemblySingleCpu {
public:
    int concurentThreadsSupported;
    std::vector<Pair> subVector;

    AssemblyParCpu(Sparse&, int);
    AssemblyParCpu(Sparse&);
    ~AssemblyParCpu();
    void simulation_per_thread();
    void eachSort(int);
    void sort() override;
    void calculateAssembly() override;
    void merge(size_t ii);
};


#endif

