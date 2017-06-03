% Parameters for Mendoza aund Yue
% If no value is known, I set to zero for now. Values come from Table II, p.916
% Default values in ().
clear all

toler = 10^(-8); %convergence parameter. Not given by MY so I set it.
bet = 0.88;
rho_e = 0.95; % autocorr of TFP shocks
sige = 0.017; % std dev of TFP shocks;
numz = 5; % (25) number of TFP realizations for Tauchen's method of TFP discretization
numb = 11; % number of bonds in bond space. Not given by MY, so I set to 11, so the 6th one is 0.
zero = 6; % index of where b = 0 in B.
cover = 3; % coverage parameter for Tauchen. This value MY don't tell, so I set to default value of 3
alpha_m = 0.43; % shares of inputs into finals' PF
alpha_l = 0.40;
alpha_k = 0.17;
lam = 0.62; % armington weight
mu = 0.65; % param of elasticity: Armington elasticity of substi is abs(1/(mu-1))
nu = 0.59; % Dixit-Stiglitz elasticity is abs(1/(nu-1)) over imported inputs j
a = 0.31;% tfp param of intermediate prod, exog constant
gam = 0.7; % CD param of intermediate prod.
thet = 0.7; % fraction of imported inputs that need to be financed using working capital
kappa = 0; % working capital loans, contracted at risk-free rate rstar
r = 0.01; % riskfree world rate, exogenous and time-invariant.
phi = 0.083; % prob. of reentering credit markets in default, exog.
omega = 1.455; % Frisch param: Frisch el = 1/(omega-1)
sigm = 2; % CRRA param
lbb = -10; % lower bound for bond space
ubb = 10; % uppe bound for bond space

k = 2; % assume capital is constant at 1. MY don't say what value they assume.
% apparently, the sol of factor markets is sensitive to this, so experiment
% around

% integration of price of imported input (assuming pstar(j) = 1 and equal
% pjstar(i) for all i,j.
pstar = (thet*(1+r)^(nu/(nu-1)) + 1 + thet)^((nu-1)/nu); % w/o default
p_aut = (1+thet)^((nu-1)/nu); % w/ default

save params