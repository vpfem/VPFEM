#include "AssemblySingleCpu.h"

AssemblySingleCpu::AssemblySingleCpu(Sparse& s) 
: Assembly(s)
{
    INFO("Assembly single CPU created");
    // HOST vectors for single CPU and par CPU
    v_value.assign(stiffMat.value, stiffMat.value + stiffMat.valueSize);
    v_dofi.assign(stiffMat.i, stiffMat.i + stiffMat.valueSize);
    v_dofj.assign(stiffMat.j, stiffMat.j + stiffMat.valueSize);
    
    // Build the host vector of indices
    h_indices.reserve(stiffMat.valueSize);
    for (size_t i = 0; i < stiffMat.valueSize; i++)
        h_indices.push_back(i);
}

AssemblySingleCpu::~AssemblySingleCpu() 
{
    INFO("Assembly single CPU destroyed");
    v_value.clear();
    v_dofi.clear();
    v_dofj.clear();
    h_indices.clear();
    h_rowPtr.clear();
}

void AssemblySingleCpu::sort()
{
    Timer timer("Time spend in sort -> ");
    std::sort(h_indices.begin(), h_indices.end(), sort_indices_assembly(v_dofi,v_dofj));
}

void AssemblySingleCpu::nnzFinder()
{
    h_rowPtr.reserve(stiffMat.get_numberOfRows()+1);
    h_rowPtr.push_back(0);
    size_t rowCounter = 1;
    size_t pointer = 0;
    for (size_t i = 0; i <= v_value.size() && v_value[h_indices[i-1]] != 0; i++)
    {
        if (rowCounter != v_dofi[h_indices[i]]) {
            h_rowPtr.push_back(pointer);
            rowCounter++;
        }
        pointer++;   
    }
}

void AssemblySingleCpu::addDuplicates() 
{
    std::vector<size_t> newZeros(h_rowPtr.size(),0);
    for (size_t row = 0; row < h_rowPtr.size(); row++) {
        newZeros[row+1] = newZeros[row] + addDuplicatesEachRow(row);
    }
    for (size_t i = 0; i < h_rowPtr.size(); i++)
    {
        h_rowPtr[i] = h_rowPtr[i] - newZeros[i];
    }
    // put zeros that accured in the system at the end of the index vector
    std::stable_partition(h_indices.begin(), h_indices.end(), [this](size_t n){return v_value[n]!=0;});
}

size_t AssemblySingleCpu::addDuplicatesEachRow(size_t row)
{
    size_t counter;
    size_t newZero_inRow = 0;
    for (size_t i = h_rowPtr[row] ; i < h_rowPtr[row+1]-1; ++i)
    {
        counter = 1;
        while (v_dofj[h_indices[i]] == v_dofj[h_indices[i+counter]] && i+counter < h_rowPtr[row+1]) {
            v_value[h_indices[i]] = v_value[h_indices[i]] + v_value[h_indices[i+counter]];
            v_value[h_indices[i+counter]] = 0;
            newZero_inRow++;
            counter++;
        }
        if (abs(v_value[h_indices[i]]) < tolorance)
        {
            v_value[h_indices[i]] = 0;
            newZero_inRow++;
        }
        i = i + counter - 1;
    }
    return newZero_inRow;
}

void AssemblySingleCpu::eraseDuplicands()
{
    v_value.erase(v_value.begin() + h_rowPtr.back() , v_value.end());
    v_dofi.erase (v_dofi.begin()  + h_rowPtr.back() , v_dofi.end());
    v_dofj.erase (v_dofj.begin()  + h_rowPtr.back() , v_dofj.end());
}


void AssemblySingleCpu::calculateAssembly()
{
    Timer *t = new Timer("Time in assembly -> ");
    sort();
    nnzFinder();
    addDuplicates();
    eraseDuplicands();
    apply_permutation_in_place(v_value, h_indices);
    apply_permutation_in_place(v_dofi, h_indices);
    apply_permutation_in_place(v_dofj, h_indices);
    delete t;
    stiffMat.valueSize = h_rowPtr.back();
    for (size_t i = 0; i < stiffMat.valueSize; i++) {
        stiffMat.value[i] = v_value[i];
        stiffMat.i[i] = v_dofi[i];
        stiffMat.j[i] = v_dofj[i];
    }
}

// ---------------------- sort_indices struct -------------------------
sort_indices_assembly::sort_indices_assembly(std::vector<size_t>& i, std::vector<size_t>& j)
  : dof_i(i), dof_j(j) {};
bool sort_indices_assembly::operator()(size_t i, size_t j) const 
{ 
    if (dof_i[i] == 0 )
        return false;
    if (dof_i[j] == 0)
        return true;
    if (dof_i[i] == dof_i[j])
        return dof_j[i] < dof_j[j];
    return dof_i[i] < dof_i[j];
};