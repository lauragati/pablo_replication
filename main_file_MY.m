% The main file for Laura Gati's replication of Mendoza & Yue (2012)
% Methods Class of Pablo Guerron, Spring Semester 2017

clear all; 
clc

% Note: This file has been kept simple to give the user a quick overview  
% of the layout. The subfiles are extensively commented to give insight
% into the details.

% Read in parameters
parameters_MY

% Discretize TFP and debt space and do a little renaming to make your life easier
[Z,P,b] = tauchen_MY(cover,sigma_z, rho_z, mu_z, numz, lbb, ubb, numb); 
e=exp(Z'); % shock process
E = size(e,2); % =numz, the size of TFP space; rename for simplicity
B = numb; % the size of debt space; rename for simplicity

% Factor market equilibrium
% As Mendoza & Yue, I first solve the static factor allocation problem and
% feed it as an input later to the VFI. To solve for the factor market
% equilibrium, I simplify the 8-dimensional system to a 3-dimensional one
% in the two types of labor and mstar, the imported input.
factor_market_equilibrium_MY

% Solve the model using VFI
vfi_MY
% This generates two figures (Fig 1 and 2) which show the relationship between bond price
% and the level of debt and TFP.

% Simulate model and recreate Table III and Figure VI
create_figs_MY
% Table III is shown in the Command Window and my reproduction of Fig. VI
% is Figure 3.