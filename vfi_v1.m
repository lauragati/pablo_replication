%% A first version of the Mendoza and Yue VFI

% Parameters
%(If no value is known, I set to zero for now. Values come from Table II, p.916)
toler = 10^(-8); %convergence parameter. Not given by MY so I set it.
bet = 0.88;
rho_e = 0.95; % autocorr of TFP shocks
sige = 0.017; % std dev of TFP shocks;
numz = 25; % number of TFP realizations for Tauchen's method of TFP discretization
numb = 10; % number of bonds in bond space. Not given by MY, so I experiment around my default of 10.
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

%% Tauchen discretization of TFP
%Preallocating matrices to speed up the code
Z = zeros(numz,1); %discrete spanning space of TFP (I'm calling TFP z)
F = zeros(numz,numz); %Markov chain transition probability matrix

lbz = - cover*sige/sqrt((1-rho_e^2)); %lower bound of the z grid
ubz = mu + cover*sige/sqrt((1-rho_e^2)); %upper blund of the z grid
dist = ubz - lbz; %overall distance to be covered by the z grid
stepsize = dist/(numz-1); %size of each bin of the z grid
for i = 1:numz
    Z(i) = lbz + (i-1)*stepsize; %linearly spaced points
end
for i = 1:numz-1
    m(i) = (Z(i) + Z(i+1))/2; %midpoints for the discrete normal CDF
end
for j = 1:numz
    a(j) = rho_e*Z(j); %this is the conditional mean to be in the previous period
    for i = 1:numz
        if i == 1
            F(i,j) = cdf('norm',m(i),a(j),sige); %first point - lowest z'
        elseif i == numz;
            F(i,j) = 1 - cdf('norm',m(i-1),a(j),sige); %last point - highest z'
        else
            F(i,j) = cdf('norm',m(i),a(j),sige) - cdf('norm',m(i-1),a(j),sige); %middle points
        end
    end
end

% Discretization for bonds (b)
lbb = -10; % lower bound for bonds (max borrowing)
ubb = 10; %upperbound for bonds (max lending) - need to check its own value in steady state to get an idea
% These two are not given by MY, so I set some symmetric bounds and will
% change values. Start with [-10,10]
distb = ubb - lbb; %overall distance to be covered
stepsize_b = distb/(numb-1); %size of each bin
for i = 1:numb
    B(i) = lbb + (i-1)*stepsize_b; %linearly spaced points in bonds space
end

%% Preparing the Big Loop to get the convergence for the value function
index = 0; %initial value to count the number of iteration
err = 1/eps; %initial value for the error. Just a pre-allocation
v = zeros(numb,numz); %initial guess for the value function. 
% % % Now I do not
% % % %need it since I have saved v from previous iteration with beta = 0.1 and
% % % %rho = 0.5.
% % % load v %This is a saved v from previous iterations where beta and rho where low. It has
% % % %a smart way to speed up the convergence
vnew = zeros(numb,1); %preallocating the updated value function
vnew_max = zeros(numb,numz); %matrix to store the maximized vnew over k'
l = zeros(numb,numz,numb); %optimal rule for agg labor given b, b', and z
c = zeros(numb,numz,numb); %optimal rule for consumption given l, b, b', and z

p = 0; % initial probability of default
%Starting the big loop to find convergence
while err > toler*max(max(abs(v))) && index < 50; %if you dont get convergence in 50 iters: stop!
    %*max(max(abs(v)))
    for b = 1:numb %loop for a given b (endogenous state variable)
        for z = 1:numz %loop for a given z (exogenous state variable)
            for bp = 1:numb %loop for a given b' (choice variable)
                
                %1.) For an initial probability p0 (to be updated in the
                %loop), get the price of borrowing q
                if B(b) >= 0 % if we're lenders
                q = 1/(1+r);
                elseif B(b) < 0 % if we're borrowers, q adjusts for risk
                q = (1-p)/(1+r);
                end
                %2.) Integrate CES to get pstar (p. 901) Calculate both
                %pstar (no default) and p_aut (default)
                %(for the default case).
                pstar = 3; % stand-in
                p_aut = 2; % 
                %3.) solve for equilibrium in factor markets; that is, given TFP
                %realization, risk-free int. rate r and pstar, get 
                % [M mstar md lm lf l pm w] (p. 903) Also do this twice,
                % for the no-default and default case.
                M     = 2; % stand-in values for now
                mstar = 3;
                md    = 4;
                lm    = 2;
                lf    = 1;
                l     = 3;
                pm    = 1;
                w     = 1;
                
                M_d     = 2; % stand-in values for now
                mstar_d = 0;
                md_d    = 5;
                lm_d    = 3;
                lf_d    = 2;
                l_d     = 5;
                pm_d    = 2;
                w_d     = 0.9;
                %4.) Use RC in eq 23 to get c (p. 906) and c_d (default)
                f_nd = Z(z)*M^alpha_m*lf^alpha_l*k^alpha_k; % prod fuct of final producer
                f_d  = Z(z)*M_d^alpha_m*lf_d^alpha_l*k^alpha_k; % PF in case of default
                c_nd(b,z,bp) = - q*B(bp) + B(b) + Z(z)*f_nd - mstar*pstar; % cons w/o default
                c_d(b,z,bp)  =Z(z)*f_d -mstar_d*p_aut; % cons w/ default
                %5.) Use gov BC to find transfer t (is the same for default) (p.908)
                t(b,z,bp) = q*B(bp) - B(b);
                %6.) Given these stuff, calculate u_nd and u_d
                u_nd(b,z) = ((c_nd(b,z,bp)-l^omega/omega)^(1-sigm)-1)/(1-sigm);
                u_d(b,z) = ((c_d(b,z,bp)-l^omega/omega)^(1-sigm)-1)/(1-sigm);
                % --> now go outside the b_p loop to calculate v_nd and v_d
                
                
                oldstuff = 0; % skip marco's stuff that is here for planning purposes
                if oldstuff == 1
                options = optimoptions('fsolve','display','off'); %option to stop fsolve to show...
                %results for each l in order to speed up the code
                foc = @(l)((Z(z)*bigk(b)^(alpha)*l^(1-alpha)...
                    + (1-delta)*bigk(b) - bigk(bp))^(-1)*...
                    (1-alpha)*Z(z)*bigk(b)^(alpha)*l^(-alpha)...
                    + psi*l^(eta-1)); % *U_c(.)F_l(.,.,.) - U_l(.) - Optimal rule for labor
                l(bp) = fsolve(foc,1/5,options); %matlab function to solve for the optimal rule for l
                l(bp) = real(l(bp)); %sometimes fsolve looks for complex number. I want to avoid it!
                l(bp) = max(l(bp),0.1); %labor should be between zero and 1. I impose nonnegativity
                c(bp) = Z(z)*(bigk(b)^(alpha))*(l(bp)^(1-alpha)) + (1-delta)*bigk(b) - bigk(bp); 
                %optimal rule for c given l, k, k', and z
                c(bp) = max(c(bp),0.1); %I want to avoid nonnegativity of consumption
                vnew(bp) = log(c(bp)) - psi*(l(bp)^eta)/eta + beta*v(b,z); %Here, I am defining
                %the dynamic programming problem for given k, k', and z +
                %optimal rule for for labor
                end
            end
            
            % evaluate 'temporary' Vnd and Vd
            v_nd_tmp(b,z) = u_nd(b,z) + bet*v(b,z); % i'm ignoring the Expectation in front of v for now.
            v_d_tmp(b,z)  = u_D(b,z) + bet*(1-phi)*v_d(tfp_v) + bet*phi*E v(0,tfp_p) %????
            v_nd_new = max(v_nd_tmp) %CHECK!!
            %7.) Determine D = 1 if v_d > v_nd, else 0
                %8.) Set vnew(bp) = max(v_d, v_nd)
                %9.) Update p as the integral over D of F (eq 27) (p. 907).
                %10.) Do the max: vnew_max(b,z) = max(vnew) and go on to
                %while loop to check error convergence
            vnew_max(b,z) = max((vnew)); % for every k and z, pick the...
            %the maximum value of vnew over k'. I will do it for each k and
            %z in the grid defined above
        end
    end
    err_matrix = vnew_max - v; %This is a matrix of the difference between v and vnew
    err = max(max(abs(err_matrix))) %The error for the convergence will be the largest...
    %in absolute value
    v = vnew_max; %Set v equal to vnew to restart the loop untill convergence
    index = index + 1; %add one value to the index to stop the loop if convergence is not reached
    %before 50 iterations.
end