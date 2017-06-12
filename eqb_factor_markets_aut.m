function [s1, s2, s3, s4, s5, s6, s7, s8] = eqb_factor_markets(z)
% % This file sets up and solves the factor market equilibrium for the
% autarky case. (Pstar gets substituted by p_aut.)

load params
syms M mstar md  lm lf l pm w

k = 1e+2; % to edit

f(1) = alpha_m*z*k^alpha_k * M^(alpha_m - mu)*lf^alpha_l*(1-lam)*mstar^(mu-1) - p_aut;
f(2) = -M + (lam*md^mu + (1-lam)*mstar^mu)^(1/mu);
f(3) = alpha_m * z * k^alpha_k * M^(alpha_m - mu) * lf^alpha_l * lam * md^(mu-1) - pm;
f(4) = alpha_l * z * k^alpha_k * M^alpha_m * lf^(alpha_l -1) - w;
f(5) = gam * pm * a * lm^(gam-1) - w;
f(6) = l^(omega-1) - w;
f(7) = lf + lm -l;
f(8) = a * lm^gam -md;

vars = [M, l, lf, lm, md, mstar, pm, w];
range = [0 5; 0 1; 0 1; 0 1; 0 1; 0 1; 0 5; 0 5];
      %   M    l    lf   lm   md mstar  pm   w
range2 = [-200, 200; -100 100; -100 100; -1 1; -1 1; -800 800; -30 30; -10 10];
      %   M             l           lf   lm     md      mstar     pm       w
range3 = [-Inf, 200; -Inf 100; -Inf 100; -Inf 1; -Inf 1; -Inf 800; -Inf 30; -Inf 10];
      %   M             l           lf   lm     md      mstar     pm       w
range4 = [50, 75;    26 34;     26 34;   0.2 0.3;     0.05 0.2;    250 310;    13 15;    4 5];
      %   M             l           lf     lm             md         mstar       pm       w      
%tic
S = vpasolve(f, vars, range4); 
% toc
 

sol_numeric = double(struct2array(S));
% The order is: [M, l, lf, lm, md, mstar, pm, w]
s1 = sol_numeric(1); % M
s2 = sol_numeric(2); % l
s3 = sol_numeric(3); % lf
s4 = sol_numeric(4); % lm
s5 = sol_numeric(5); % md
s6 = sol_numeric(6); % mstar
s7 = sol_numeric(7); % pm
s8 = sol_numeric(8); % w
