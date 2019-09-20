#include <iostream>
#include "Application/Application.h"

int main(int argc, char **argv) 
{
    Timer timer;
    Application app(atol(argv[1]), atol(argv[2]));
    app.ModelBuilder();
    app.MaterialBuilder();
    timer.Start();
    app.ElasticRun();
    timer.Stop();

    double* d = app.get_displacement();
    int nx = atol(argv[1]);
    int ny = atol(argv[2]);
    //std::cout << "Bottom Disp: "  << d[nx*2-1] << ", Top disp: " << d[(ny+1)*nx*2-1] << std::endl;
    //std::cout << d[nx*2-1] << std::endl;

}