#include "Assembly.h"

Assembly::Assembly(Sparse& s)
: stiffMat(s)
{
    INFO("Assembly created");
    tolorance = 0.001;

    // DEVICE vectors for GPU
    //d_value.assign(stiffMat.value, stiffMat.value + stiffMat.valueSize);
    //d_dofi.assign(stiffMat.i, stiffMat.i + stiffMat.valueSize);
    //d_dofj.assign(stiffMat.j, stiffMat.j + stiffMat.valueSize);



    // for (auto& ii : h_indices) {
    //     std::cout<< ii << "\n";
    // }
}

Assembly::~Assembly() 
{
    INFO("Assembly distroyed");
}


// -- override the cout << oprator 
std::ostream& operator<< (std::ostream &out, Assembly const& assemble) 
{
    for (size_t i = 0; i < assemble.stiffMat.valueSize; i++) {
        out << '\t' << assemble.stiffMat.i[i] << '\t' << assemble.stiffMat.j[i] << '\t' << assemble.stiffMat.value[i] << "\n";
    }
    return out;
}