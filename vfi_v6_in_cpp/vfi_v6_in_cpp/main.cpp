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



/* Includes, system */
#include <iomanip>
#include <fstream>
#include <cmath> // apparently needed for sqrt, fab, exp. Pablo didn't need it.


using namespace std;

// A.) DEFINE SOME PARAMETERS (those which are needed as double in code body are commented out here, but kept here for readability)
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

// B.) DEFINE TAUCHEN'S DISCRETIZATION METHOD


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
    
    // (3.) FACTOR MARKET EQUILIBRIUM
    
    // (4.) VFI
    
    // (5.) SIMULATION
    // recreate HP-filter using Eigen library here.
    
    /* Next up:
     [Z,P,b] = tauchen_MY(cover,sigma_z, rho_z, mu_z, numz, lbb, ubb, numb);
     e=exp(Z'); % shock process
     E = size(e,2); % =numz (rename for simplicity)
     */
    
}
