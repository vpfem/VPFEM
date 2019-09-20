#include "ElasticPlainStrain.h"

ElasticPlainStrain::ElasticPlainStrain(double E, double nu)
: ElasticModulus(E), poissonRatio(nu)
{
    INFO("Material created");
}

void ElasticPlainStrain::elasticMaterialMatrix(double* materialMatrix, double* secantModulus) 
{
    // this function will create a vector for material matrix
    // D11 D22 D33 D12 D13 D23
    double E1 = ElasticModulus*(1-poissonRatio)/((1-poissonRatio*2)*(1+poissonRatio));
    double E2 = poissonRatio*E1/(1-poissonRatio);
    double G  = ElasticModulus/(2*(1+poissonRatio));
    materialMatrix[0] = E1;
    materialMatrix[1] = E1;
    materialMatrix[2] =  G;
    materialMatrix[3] = E2;
    materialMatrix[4] = 0;
    materialMatrix[5] = 0;

    // store the secant stiffness
    secantModulus[0] = E1; //D11 which is Ec1
    secantModulus[1] = E1; //D22 which is Ec2
    secantModulus[2] = 0.0; //Ds_x which is Es1
    secantModulus[3] = 0.0; //Ds_y which is Es2
}

ElasticPlainStrain::~ElasticPlainStrain()
{
}