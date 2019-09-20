#ifndef ASSEMBLYSINGLECPU_H
#define ASSEMBLYSINGLECPU_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <thrust/sort.h>

#include "../Log/Log.h"
#include "../Timer/Timer.h"
#include "../Sparse/Sparse.h"
#include "Assembly.h"

class AssemblySingleCpu : public Assembly {
public:
    std::vector<double> v_value;
    std::vector<size_t> v_dofi;
    std::vector<size_t> v_dofj;

    std::vector<size_t> h_indices;

    std::vector<size_t> h_rowPtr;
private:

public:
    AssemblySingleCpu(Sparse& s);
    ~AssemblySingleCpu();
    
    void sort() override;
    void calculateAssembly() override;
    void nnzFinder();
    void addDuplicates();
    size_t addDuplicatesEachRow(size_t) ;
    void eraseDuplicands();


    template <typename T>
    void apply_permutation_in_place(
    std::vector<T>& vec,
    const std::vector<std::size_t>& p)
    {
        std::vector<bool> done(vec.size());
        for (std::size_t i = 0; i < vec.size(); ++i)
        {
            if (done[i])
            {
                continue;
            }
            done[i] = true;
            std::size_t prev_j = i;
            std::size_t j = p[i];
            while (i != j)
            {
                std::swap(vec[prev_j], vec[j]);
                done[j] = true;
                prev_j = j;
                j = p[j];
            }
        }
    }

private:

};


// ---------------------- sort_indices struct -------------------------
struct sort_indices_assembly {
private:
    std::vector<size_t>& dof_i; // the variable that index is going to sorted base of
    std::vector<size_t>& dof_j;
public:
    sort_indices_assembly(std::vector<size_t>&, std::vector<size_t>&);
    bool operator()(size_t i, size_t j) const;
};
#endif