% This file is more a note of Mendoza & Yue than a working file. It
% outlines the parameters, their values and thus tries to form an idea of
% what the VFI should look like.

% Agents
% 1. domestic HHs: 
% choose C and L, do not borrow or lend, take wages, profits
% and gov transfers as given
% 2. domestic firms
% - final goods producers (f sector)
% final good = CD of TFP, domestic intermediates (md), imported intermediates (mstar) and
% labor in f (lf).
% choose kappa so that eq. 6, the collateral constraint, holds with
% equality
% - intermediate producers (md sector)
% produce md using labor (lm)
% 3. domestic gov
% the sovereign chooses whether to default or not and chooses the level of
% borrowing, b(t+1),
% cannot commit to repayment
% when defaults, excluded from credit markets right away and reenters with
% prob. theta. NEW: in this case, also final producers are excluded :D
% 4. foreign lenders
% provide working capital loans, and I guess buy the gov bond

% --> whole thing is set up as a planner problem where the gov chooses
% debt policy (b and whether to default)
% consumption
% factor allocations
% States are (b tfp), and taken as given is q(b_p,tfp)

% Variables (If no value is known, I set to zero for now. Values come from Table II, p.916)
bet = 0.88;
rho_e = 0.95; % autocorr of TFP shocks
sige = 0.017; % std dev of TFP shocks;
n = 25; % number of TFP realizations for Tauchen's method of TFP discretization
f = tfp*M^alpha_m*lf^alpha_l*k^alpha_k; % prod fuct of final producer
alpha_m = 0.43; % shares of inputs into finals' PF
alpha_l = 0.40;
alpha_k = 0.17;
M = (lam*md^mu + (1-lam)*mstar^mu)^(1/mu); % Armington aggregator combining domestic inputs and imported ones
lam = 0.62; % armington weight
mu = 0.65; % param of elasticity: Armington elasticity of substi is abs(1/(mu-1))
mstar = 0; % imported input aka imported intermediate, NEED TO TAKE INTEGRAL over imported inputs j
nu = 0.59; % Dixit-Stiglitz elasticity is abs(1/(nu-1)) over imported inputs j
a = 0.31;% tfp param of intermediate prod, exog constant
gam = 0.7; % CD param of intermediate prod.
md = 0; %domestic input aka domestic intermediate
pm = 0; % relative price of domestic input to imported input, endogenous
pjstar = 0; %price of imported input j, constant and exogenous.
pstar = 0; % price of agg imported input, the CES of imported inputs of price pjstar and 
% those imported inputs that need financing, so their price is augmented by rstar
p_aut = 0; % price of agg imported inputs in autarky, the CES of only non-financed imported inputs
thet = 0.7; % fraction of imported inputs that need to be financed using working capital
kappa = 0; % working capital loans, contracted at risk-free rate rstar
r = 0.01; % riskfree world rate, exogenous and time-invariant.
b = 0; % borrowing = gov bonds, the ammount to be paid tomorrow. <0 for borrowing. b belongs to [bmin bmax]
bmin = 0;
bmax = 0;
phi = 0.083; % prob. of reentering credit markets in default, exog.
q = 0; %price of bond tomorrow; q(b_p,tfp)
v = 0; %overall value function of the gov
v_d = 0; %value in case of default
c = 0; % consumption
l = 0; % labor
lf = 0; %labor in final prod
lm = 0; % labor in domestic prod of intermediate
omega = 1.455; % Frisch param: Frisch el = 1/(omega-1)
sigm = 2; % CRRA param
u = ((c-l^omega/omega)^(1-sigm)-1)/(1-sigm);
v_nd = max(u + bet*Ev_p); %value in case of no default % NEED TO TAKE E
c + q*b_p - b = tfp*f - mstar*pstar; % st 1
l = lf + lm; % st 2
md = a*lm^gam; % st 3
v_d = max(u + bet*(1-phi)E v_d(tfp_v) + bet*phi*E v(0,tfp_p)); % NEED TO TAKE E
c = tfp*f -mstar*p_aut % st 1. (new)
l = lf + lm; % st 2 (same as in no default)
md = a*lm^gam; % st 3 (same as in no default)

v = max(v_d, v_nd);

D = 0; %default set. =1 if default, zero otherwise.
z = 0; % for each value of tfp, markov transition probability of tfp_p.
p = 0; % prob of default. NEED to take integral of z for all tfp such that we default (D=1).

q = 1/(1+r) % if b_p > 0
q = 1-p/(1+r) % if b_p <0

% the VFI will have to:
% 1.) discretize the debt space
% 2.) discretize the TFP space
% 3.) for each inherited debt position, for each expected TFP level, determine the probability of default by
% comparing the gov's value when defaulting vs. when staying alive