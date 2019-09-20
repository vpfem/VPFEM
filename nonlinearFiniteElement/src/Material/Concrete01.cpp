#include "Concrete01.h"

#define PI 3.14159265359

Concrete01::Concrete01(double fcPrime, double fcr, double Ec, double e0, 
                        double nu, double Es, double beta,
                        double ro_x, double fyx, double alpha_x, 
                        double ro_y, double fyy, double alpha_y)
:  ElasticPlainStress(Ec, nu), fcPrime(fcPrime), fcr(fcr), Ec(Ec), e0(e0), 
    nu(nu), Es(Es), beta(beta), ro_x(ro_x), fyx(fyx), alpha_x(alpha_x), 
    ro_y(ro_y), fyy(fyy), alpha_y(alpha_y)
{
    e_cr = fcr*e0/(2*fcPrime);
}

Concrete01::~Concrete01()
{

}
double Concrete01::errorCheck(double E1, double E2, double Es1, double Es2, double* secantModulus)
{
    double _error[4];
    _error[0] = std::abs((secantModulus[0] -  E1) / E1 );
    _error[1] = std::abs((secantModulus[1] -  E2) / E2 );
    _error[2] = std::abs((secantModulus[2] - Es1) / Es1);
    _error[3] = std::abs((secantModulus[3] - Es2) / Es2);
    secantModulus[0] = E1;
    secantModulus[1] = E2;
    secantModulus[2] = Es1;
    secantModulus[3] = Es2;
    return *std::max_element(_error, _error+4);
}
// This function calculate material matrix using new strain
double Concrete01::UpdateMatMatrix(double* materialMatrix, double* secantModulus, double* e) 
{
    // Find principal strains
    double ec[2]; double theta;
    PrincipalStrain(ec, theta, e);
    double es[2];
    es[0] = SteelStrain(e, beta, alpha_x);
    es[1] = SteelStrain(e, beta, alpha_y);
    // Determine the principal stresses
    double fc[2]; StrainToStressConcerete(fc, ec);
    double fs[2]; StrainToStressSteel(fs, es); // 0-x 1-y
    // Determine the secant stiffness
    // 0-D11 1-D22 2-D33 3-D12 4-D13 5-D23
    double Dc[6] = {0}; matMatrixConcrete(Dc, ec, fc);
    double Ds_x[6] = {0}; double Ds_y[6] = {0};
    matMatrixSteel(Ds_x, Ds_y, es, fs); // 0-x 1-y
    // check error for Ec1, Ec2, Esx, Esy 
    double error = errorCheck(Dc[0], Dc[1], Ds_x[0], Ds_y[0], secantModulus);
    // Transform the material matrices
    transformation(Dc, theta);
    if (ro_x != 0) {transformation(Ds_x, alpha_x + beta);}
    if (ro_y != 0) {transformation(Ds_y, alpha_y + beta);}
    // Add secant stiffnesses together
    for (int i = 0; i < 6; i++) { materialMatrix[i] = Dc[i] + Ds_x[i] + Ds_y[i];}
    return error;
}

void Concrete01::PrincipalStrain(double* result, double& theta, double* strain)
{
    theta = PI/2 + 0.5 * atan (strain[2]/(strain[0]-strain[1])); // The reason to adding PI/2 is that we are looking at the angle with a diffrent axis. 
    result[0] = (strain[0] + strain[1]) / 2 + 0.5 * pow((((strain[0]-strain[1]) * (strain[0]-strain[1])) + (strain[2]*strain[2])) , 0.5);
    result[1] = (strain[0] + strain[1]) / 2 - 0.5 * pow((((strain[0]-strain[1]) * (strain[0]-strain[1])) + (strain[2]*strain[2])) , 0.5);
}

double Concrete01::SteelStrain(double* strain, double beta, double alpha)
/*
beta: angle of the element
alpha: angle of the steel rebars
*/
{
    double theta = beta + alpha;
    return (strain[0]+strain[1])/2 + ((strain[0]-strain[1])/2)*cos(2*theta) + sin(2*theta)*strain[2]/2;
}

void Concrete01::StrainToStressConcerete(double *result, double* strain_p) 
{
    if (strain_p[0] <= e_cr) 
    {
        result[0] = 2*fcPrime*strain_p[0]/e0;
    } else if (strain_p[0] > e_cr)
    {
        result[0] = fcr/(1+sqrt(200*strain_p[0]));
    }
    
    double fc2max = std::max(-fcPrime /(0.8 - (0.34*strain_p[0]/e0)),-fcPrime);
    result[1] = fc2max * (2.0*(strain_p[1]/e0) - pow((strain_p[1]/e0),2));
}

void Concrete01::StrainToStressSteel(double* result, double* strain_p) 
{   
    if (strain_p[0] <= 0)
    {
        result[0] = std::min(Es*strain_p[0],fyx);
    }
    else
    {
        result[0] = std::max(Es*strain_p[0],-fyx);
    }
    // for steel in y direction
    if (strain_p[1] <= 0)
    {
        result[1] = std::min(Es*strain_p[1],fyy);
    }
    else
    {
        result[1] = std::max(Es*strain_p[1],-fyy);
    }
}

void Concrete01::matMatrixConcrete(double* mat, double* strain, double* stress)
{
    // 0-D11 1-D22 2-D33 3-D12 4-D13 5-D23
    mat[1] = stress[0]/strain[0];
    mat[0] = stress[1]/strain[1];
    mat[2] = mat[0] * mat[1] / (mat[0] + mat[1]);
}

void Concrete01::matMatrixSteel(double* mat_x, double* mat_y, double* strain, double* stress)
{
    mat_x[0] = ro_x * stress[0] / strain[0];
    mat_y[0] = ro_y * stress[1] / strain[1];
}

void Concrete01::transformation(double* mat, double angle)
{
    double c = std::cos(angle);
    double s = std::sin(angle);
    double c4 = std::pow(c,4);
    double s4 = std::pow(s,4);
    double c3 = std::pow(c,3);
    double s3 = std::pow(s,3);
    double c2 = std::pow(c,2);
    double s2 = std::pow(s,2);

    double result[6];
    result[0] = mat[0]*c4 + 2*mat[3]*c2*s2 - 4*mat[4]*c3*s + mat[1]*s4 - 4*mat[5]*c*s3 + 4*mat[2]*c2*s2;
    result[1] = mat[0]*s4 + 2*mat[3]*c2*s2 + 4*mat[4]*c*s3 + mat[1]*c4 + 4*mat[5]*c3*s + 4*mat[2]*c2*s2;
    result[2] = c*s*(mat[0]*c*s - mat[3]*c*s + mat[4]*(c2 - s2)) - c*s*(mat[3]*c*s - mat[1]*c*s + mat[5]*(c2 - s2)) + (c2 - s2)*(mat[4]*c*s - mat[5]*c*s + mat[2]*(c2 - s2));
    result[3] = c2*(mat[3]*c2 + mat[1]*s2 - 2*mat[5]*c*s) + 2*c*s*(mat[4]*c2 + mat[5]*s2 - 2*mat[2]*c*s) + s2*(mat[0]*c2 + mat[3]*s2 - 2*mat[4]*c*s);
    result[4] = c*s*(mat[0]*c2 + mat[3]*s2 - 2*mat[4]*c*s) - c*s*(mat[3]*c2 + mat[1]*s2 - 2*mat[5]*c*s) + (c2 - s2)*(mat[4]*c2 + mat[5]*s2 - 2*mat[2]*c*s);
    result[5] = c*s*(mat[0]*s2 + mat[3]*c2 + 2*mat[4]*c*s) - c*s*(mat[3]*s2 + mat[1]*c2 + 2*mat[5]*c*s) + (c2 - s2)*(mat[4]*s2 + mat[5]*c2 + 2*mat[2]*c*s);

    std::copy(result, result+6, mat);
}