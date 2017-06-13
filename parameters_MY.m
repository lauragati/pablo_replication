% Parameters for Mendoza & Yue (2012)

% Unless otherwise indicated, all values are taken from MY, baseline
% calibration, Table II, p. 916.
alpha_m=.43; % share of intermediate goods in gross output
gam=.7; % labor income share in value added and in intermediate goods production
alpha_L=(1-alpha_m)*gam;
alpha_k=1-alpha_m-alpha_L;
sigm=2; % coef. of relative risk aversion
r=1.01; % risk-free interest rate
omega=1.455; % curvature of labor supply, 1/(omega-1)=Frisch wage elasticity
phi=.083; % reentry probability
lam=.62; % Armington weight of domestic inputs
mu=.65; % Armington curvature parameter
nu=.59; % Dixit-Stiglitz curvature parameter
A=.31; % intermediate goods TFP coefficient
bet=.88; % subjective discount factor
thet=.7; % upper bound of imported inputs with working capital
xi=0;%-.67; % TFP semi-elasticity of exogenous capital flows
% TFP process parameters
rho_z=.95; % AR-coefficient
sigma_z=0.017; % variance
mu_z=0;%conditional mean 

numz = 25; % number of TFP realizations for Tauchen's method of TFP discretization
numb = 100; % number of bonds in bond space. Not given by MY, so I set to 100.
cover = 3; % coverage parameter for Tauchen. This value MY don't tell, so I set to default value of 3
lbb = 0; % lower bound for bond space (also not given by MY)
ubb = .05; % upper bound for bond space (also not given by MY)