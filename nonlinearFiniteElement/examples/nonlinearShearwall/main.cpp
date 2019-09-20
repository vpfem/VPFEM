#include <iostream>
#include "Application/Application.h"

int main(int argc, char **argv) 
{
    Timer timer;
    timer.Start();
    Application app(atol(argv[1]), atol(argv[2]));
    app.ModelBuilder();
    app.MaterialBuilder();
    app.NonlinearRun();
    timer.Stop();
    //app.ElasticRun();
    //ARRAY(app.get_d_matrix(),app.get_numberOFElements()*6);
    double strain[3] = {0};
    double* app_strain = app.get_strain();
    for (size_t i = 0; i < app.get_numberOFElements(); i++)
    {
        strain[0] = strain[0] + (app_strain[3*i+0])/app.get_numberOFElements();
        strain[1] = strain[1] + (app_strain[3*i+1])/app.get_numberOFElements();
        strain[2] = strain[2] + (app_strain[3*i+2])/app.get_numberOFElements();
    }
    ARRAY(strain,3);
    //ARRAY(app_strain,app.get_numberOFElements()*3);
    //ARRAY(app.get_displacement(),app.get_Geometry()->get_x_y_size()*2-3);
    //ARRAY(app.get_Geometry()->get_Load().loadVec,app.get_Geometry()->get_x_y_size()*2-3);
    //ARRAY(app.get_Geometry()->get_Dof().get_free(),app.get_Geometry()->get_x_y_size()*2);
}