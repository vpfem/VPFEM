#include "AssemblyParCpu.h"

AssemblyParCpu::AssemblyParCpu(Sparse& s) 
: AssemblyParCpu(s, std::thread::hardware_concurrency()-1)
{}

AssemblyParCpu::AssemblyParCpu(Sparse& s, int numberOfCores) 
: AssemblySingleCpu(s), concurentThreadsSupported(numberOfCores)
{
    INFO("Assembly par CPU created with " + std::to_string(concurentThreadsSupported) + " threads");
    simulation_per_thread();
}

// SubVecSize: number of subvector required
int SubVecSize(int n)
{
    if (n > 1)
        return (std::ceil(static_cast<float>(n/2.0)) + SubVecSize(n/2));
    else
        return 0;
}
// BuildSubVector: Make the subvector ready for merge
void BuildSubVector(std::vector<Pair> & s, size_t CTS)
// s is the subvector
// CTS the concurentThreadsSupported
{
    // Calculating the subvector entries
    size_t subVectorSize = CTS + SubVecSize(CTS);
    size_t counter = CTS;
    size_t max = CTS;
    size_t min = 1;
    while (counter < subVectorSize) {
        for (size_t i = min; i < max; i = i+2)
        {
            s[counter++] = Pair(std::min(s[i-1].x, s[i].x),
                                std::max(s[i-1].y, s[i].y));
        }
        min = max + 1;
        if (max%2) {
            std::rotate(s.begin() + max - 1 , s.begin() + max  , s.begin() + counter);
            min = max;
        }
        max = counter;        
    }
}

void AssemblyParCpu::simulation_per_thread()
{
    int simulationSize = v_dofi.size();
    size_t subVectorSize = concurentThreadsSupported + SubVecSize(concurentThreadsSupported);
    subVector.resize(subVectorSize);
    int loadPerThread =  static_cast<int>(simulationSize/concurentThreadsSupported);
    int loadLastThread = static_cast<int>(simulationSize/concurentThreadsSupported + simulationSize%concurentThreadsSupported);
    for (int i = 0; i < concurentThreadsSupported - 1; i++)
    {
        subVector[i] = Pair(i*loadPerThread,(i+1)*loadPerThread);
    }
    subVector[concurentThreadsSupported-1] = Pair((concurentThreadsSupported-1)*loadPerThread,(concurentThreadsSupported-1)*loadPerThread + loadLastThread);
}

AssemblyParCpu::~AssemblyParCpu() 
{
    INFO("Assembly par CPU destroyed");
}

void AssemblyParCpu::eachSort(int ii)
{
    std::sort(h_indices.begin() + subVector[ii].x, 
              h_indices.begin() + subVector[ii].y, 
              sort_indices_assembly(v_dofi,v_dofj));
}

// fix the merge to run in parallel
void AssemblyParCpu::merge(size_t ii)
{
    
    const auto comparison = [this](size_t i, size_t j) {
    if (v_dofi[i] == 0 )
        return false;
    if (v_dofi[j] == 0)
        return true;
    if (v_dofi[i] == v_dofi[j])
        return v_dofj[i] < v_dofj[j];
    return v_dofi[i] < v_dofi[j];
    };
    /*
    std::merge(h_indices.begin() + subVector[ii].x, h_indices.begin() + subVector[ii].y,
               h_indices.begin() + subVector[ii+1].x, h_indices.begin() + subVector[ii+1].y, 
               h_indices_new.begin() + std::min(subVector[ii].x, subVector[ii+1].x), 
               comparison);
    */
    std::inplace_merge(h_indices.begin() + subVector[ii].x, h_indices.begin() + subVector[ii].y,
               h_indices.begin() + subVector[ii+1].y, comparison);

    
}

void AssemblyParCpu::sort()
{
    Timer timer("Time spend in ParSort -> ");
    std::thread t[concurentThreadsSupported]; 
    for (int i = 0; i < concurentThreadsSupported; i++)
    {
        t[i] = std::thread(&AssemblyParCpu::eachSort,this,i);
    }
    for (int i = 0; i < concurentThreadsSupported; i++)
    {
        t[i].join();
    }
    BuildSubVector(subVector, concurentThreadsSupported);
    /*    for (auto& ii : subVector)
    {
        std::cout << ii << "\n";
    } */
    /*
    // Calculating the subvector entries
    size_t counter = 0;
    size_t max = concurentThreadsSupported;
    size_t min = 1;
    int it = 0;
    while (counter < subVector.size()-1) {
        it = 0;
        for (size_t i = min; i < max && counter < subVector.size()-2; i = i+2)
        {
            //std::cout<<counter<<std::endl;
            t[it] = std::thread(&AssemblyParCpu::merge,this,counter);
            counter = counter+2;
            it++;
        }
        for (int iit = 0; iit < it; iit++)
        {
            t[iit].join();
        }
        min = max;
        max = max + it;        
    }
    */


}
void AssemblyParCpu::calculateAssembly()
{
    Timer *t = new Timer("Time in assembly par cpu -> ");
    sort();
    //nnzFinder();
    //addDuplicates();
    //eraseDuplicands();
    apply_permutation_in_place(v_value, h_indices);
    apply_permutation_in_place(v_dofi, h_indices);
    apply_permutation_in_place(v_dofj, h_indices);
    delete t;
    //INFO(h_rowPtr.back());
    //stiffMat.valueSize = h_rowPtr.back();
    for (size_t i = 0; i < stiffMat.valueSize; i++) {
        stiffMat.value[i] = v_value[i];
        stiffMat.i[i] = v_dofi[i];
        stiffMat.j[i] = v_dofj[i];
    }
}
// ========================================================================
// pair class
Pair::Pair(int a, int b)
: x(a), y(b)
{}

Pair::Pair()
{
    x = 0;
    y = 0;
}

// -- override the cout << oprator 
std::ostream& operator<< (std::ostream &out, Pair const& a) 
{
    return out << "<" << a.x << ", " << a.y << ">";
}