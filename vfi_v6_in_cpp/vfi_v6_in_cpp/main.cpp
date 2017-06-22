//
//  main.cpp
//  vfi_v6_in_cpp
//
//  Created by Laura Gati on 6/22/17.
//  Copyright © 2017 Laura Gati. All rights reserved.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <tuple>
#include <iomanip>
#include <fstream>
#include <cmath> 
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

using namespace Eigen;
using namespace std;

// A.) DEFINE SOME PARAMETERS
// (those which are needed as double in code body are commented out here, but kept here for readability)
// ------
//#define alpha_m 0.43 // share of intermediate goods in gross output
//#define gam 0.7 // labor income share in value added and in intermediate goods production
//#define alpha_L 0.3990
//#define alpha_k 0.1710
//#define sigm 2 // coefficient of relative risk aversion
//#define r 1.01; // risk-free interest rate
//#define omega 1.455; // curvature of labor supply, 1/(omega-1)=Frisch wage elasticity
//#define phi 0.083; // % reentry probability
//#define lam 0.62; // Armington weight of domestic inputs
//#define mu 0.65; // Armington curvature parameter
//#define nu 0.59; // Dixit-Stiglitz curvature parameter
//#define A 0.31; // intermediate goods TFP coefficient
//#define bet 0.88; // subjective discount factor
//#define thet 0.7; // upper bound of imported inputs with working capital
//#define xi 0; //-0.67; % TFP semi-elasticity of exogenous capital flows
#define rho_z 0.95 // AR-coefficient of TFP process
#define sigma_z 0.017 // std. dev of TFP process
#define mu_z 0 // conditional mean of TFP
#define numz 25 // number of TFP realizations for Tauchen's method of TFP discretization
#define numb 100 // number of bonds in bond space. Not given by MY, so I set to 11, so the 6th one is 0.
#define cover 3 // coverage parameter for Tauchen. This value MY don't tell, so I set to default value of 3
#define lbb 0 // lower bound for bond space
#define ubb 0.05 // upper bound for bond space
//#define smooth 1600 // smoothing parameter for HP filter

// ------
// B.) DEFINE TAUCHEN'S DISCRETIZATION METHOD
// ------
// Define a Normal CDF in preamble
// A stand alone normcdf
double mynormcdf(double x) {
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
    
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    
    return 0.5*(1.0 + sign*y);
};

// Tauchen's method itself
void tauchen(double rrho, double ssigma, vector<double>& Z, vector<double>& P) {
    double ssigma_z = sqrt(pow(ssigma,2)/(1-pow(rrho,2)) );
    int nzgrid = Z.size();
    Z[nzgrid-1] = 5*ssigma_z; Z[0] = -5*ssigma_z;
    double step = (Z[nzgrid-1] - Z[0])/ double(nzgrid-1);
    for (int i = 2; i <= nzgrid-1; i++) {
        Z[i-1] = Z[i-2] + step;
    };
    
    for (int i_z = 0; i_z < nzgrid; ++i_z) {
        P[i_z] = mynormcdf( (Z[0]-rrho*Z[i_z]+step/2)/ssigma  );
        P[i_z+nzgrid*(nzgrid-1)] = 1.0 - mynormcdf( (Z[nzgrid-1]-rrho*Z[i_z]-step/2)/ssigma  );
    };
    
    for (int i_z = 0; i_z < nzgrid; ++i_z) {
        for (int i_zplus = 1; i_zplus < nzgrid-1; ++i_zplus) {
            P[i_z+nzgrid*i_zplus] = mynormcdf( (Z[i_zplus]-rrho*Z[i_z]+step/2)/ssigma  )-mynormcdf( (Z[i_zplus]-rrho*Z[i_z]-step/2)/ssigma  );
        };
    };
    
}; // Tauchen ends.

// ------
// C.) DEFINE FUNCTORS AS INPUTS TO SOLVER FOR FACTOR MARKET EQUILIBRIUM
// C.1) Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
    
    int m_inputs, m_values;
    
    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}
    
    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
    
}; // closes struct Functor

// C.2) Functor for good scenario factor market solver
struct my_functor_good : Functor<double>
{
    my_functor_good(void): Functor<double>(3,3) {}
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
    {
        // Equations 14, 15 and the one describing ms come here
        fvec(0) = 10.0*pow(x(0)+3.0,2) +  pow(x(1)-5.0,2);
        fvec(1) = 0;
        fvec(2) = 0;
        
        return 0;
    } // closes operator
}; // closes struct my_functor_good

// C.2) Functor for bad (default) scenario factor market solver




// --------------------- MAIN CODE ---------------------------------
int main(int argc, const char * argv[]) {

    // (1.) DEFINE PARAMETERS NEEDED AS DOUBLE IN BODY OF CODE
    double alpha_m = 0.43; // share of intermediate goods in gross output
    double gam = 0.7; // labor income share in value added and in intermediate goods production
    double alpha_L = 0.3990;
    double alpha_k = 0.1710;
    double sigm = 2; // coefficient of relative risk aversion
    double r = 1.01; // risk-free interest rate
    double omega = 1.455; // curvature of labor supply, 1/(omega-1)=Frisch wage elasticity
    double phi = 0.083; // % reentry probability
    double lam = 0.62; // Armington weight of domestic inputs
    double mu = 0.65; // Armington curvature parameter
    double nu = 0.59; // Dixit-Stiglitz curvature parameter
    double A = 0.31; // intermediate goods TFP coefficient
    double bet = 0.88; // subjective discount factor
    double thet = 0.7; // upper bound of imported inputs with working capital
    double xi = 0; //-0.67; % TFP semi-elasticity of exogenous capital flows
    
    // (2) TAUCHEN DISCRETIZATION
    vector<double> B(numb); // bond grid
    vector<double> Z(numz); // endowment grid
    vector<double> P(numz*numz, 0.0); // prob of default, comes from tauchen.
    
    tauchen(rho_z, sigma_z, Z, P);
    
    vector<double> e(numz, 0.0); // = exp(z), shock process in levels
    for (int i=0; i<=numz; i++){
        e[i] = exp(Z[i]);
    }
    
    // (3.) FACTOR MARKET EQUILIBRIUM
    // (3.1) Good (non-default) factor market
    Eigen::VectorXd x(3);
    x(0) = 0.07;
    x(1) = 0.05;
    x(2) = 0.005;
    std::cout << "x: " << x << std::endl;
    
    my_functor_good functor;
    Eigen::NumericalDiff<my_functor_good> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<my_functor_good>,double> lm(numDiff);
    lm.parameters.maxfev = 2000;
    lm.parameters.xtol = 1.0e-10;
    std::cout << lm.parameters.maxfev << std::endl;
    
    int ret = lm.minimize(x);
    std::cout << lm.iter << std::endl;
    std::cout << ret << std::endl;
    
    std::cout << "x that minimizes the function: " << x << std::endl;
    
    // (4.) VFI
    
    // (5.) SIMULATION
    // recreate HP-filter using Eigen library here.
    
    /* Next up:
     [Z,P,b] = tauchen_MY(cover,sigma_z, rho_z, mu_z, numz, lbb, ubb, numb);
     e=exp(Z'); % shock process
     E = size(e,2); % =numz (rename for simplicity)
     */
    
}
