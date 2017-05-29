%Closing and clearing previous projects
clear all
close all

%Calibration
beta = 0.995; %Discount factor - I set it arbitrarily high.
alpha = 1/3; %Cobb-Douglas Production Function with parameter alpha (capital share)
delta = 0.01; %Capital depreciation - arbitrarily small
rho = 0.99; %persistence of the shock to productivity
sigma = 1; %variance of the normal shock to productivity
eta = 4; % Frisch elasticity. I am assuming additive separability of the utility function
mu = 0; %growth rate of the TFP in steady state. Set it equal to zero to avoid growth in s.s.
kappa = 3; %coverage parameter from Touchen's discretization of z
numz = 10; %number of states for the discretization of z
numk = 10; %number of states for the discretization of k
toler = 10^(-8); %convergence parameter. Maybe could be bigger but I guess it is enough

%Preallocating matrices to speed up the code
bigz = zeros(numz,1); %discrete spanning space of z
F = zeros(numz,numz); %Markov chain transition probability matrix

%Touchen Discretization for z
lbz = mu - kappa*sigma/sqrt((1-rho^2)); %lower bound of the z grid
ubz = mu + kappa*sigma/sqrt((1-rho^2)); %upper blund of the z grid
dist = ubz - lbz; %overall distance to be covered by the z grid
stepsize = dist/(numz-1); %size of each bin of the z grid
for i = 1:numz
    bigz(i) = lbz + (i-1)*stepsize; %linearly spaced points
end
for i = 1:numz-1
    m(i) = (bigz(i) + bigz(i+1))/2; %midpoints for the discrete normal CDF
end
for j = 1:numz
    a(j) = (1-rho)*mu + rho*bigz(j); %this is the conditional mean to be in the previous period
    for i = 1:numz
        if i == 1
            F(i,j) = cdf('norm',m(i),a(j),sigma); %first point - lowest z'
        elseif i == numz;
            F(i,j) = 1 - cdf('norm',m(i-1),a(j),sigma); %last point - highest z'
        else
            F(i,j) = cdf('norm',m(i),a(j),sigma) - cdf('norm',m(i-1),a(j),sigma); %middle points
        end
    end
end

%Touchen Discretization for capital (k)
lbk = 0; %capital cannot be negative
ubk = 10; %upperbound for capital - need to check its own value in steady state to get an idea
distk = ubk - lbk; %overall distance to be covered
stepsizek = distk/(numk-1); %size of each bin
for i = 1:numk
    bigk(i) = lbk + (i-1)*stepsizek; %linearly spaced points
end

%Preparing the Big Loop to get the convergence for the value function
index = 0; %initial value to count the number of iteration
err = 1/eps; %initial value for the error. Just a pre-allocation
%v = zeros(numk,numz); %initial guess for the value function. Now I do not
%need it since I have saved v from previous iteration with beta = 0.1 and
%rho = 0.5.
load v %This is a saved v from previous iterations where beta and rho where low. It has
%a smart way to speed up the convergence
vnew = zeros(numk,1); %preallocating the updated value function
vnew_max = zeros(numk,numz); %matrix to store the maximized vnew over k'
l = zeros(numk,1); %optimal rule for labor given k, k', and z
c = zeros(numk,1); %optimal rule for onsumption given l, k, k', and z

%Starting the big loop to find convergence
while err > toler*max(max(abs(v))) && index < 50; %if you dont get convergence in 50 iters: stop!
    %*max(max(abs(v)))
    for i = 1:numk %loop for a given k (endogenous state variable)
        for j = 1:numz %loop for a given z (exogenous state variable)
            for k = 1:numk %loop for a given k' (choice variable)
                options = optimoptions('fsolve','display','off'); %option to stop fsolve to show...
                %results for each l in order to speed up the code
                foc = @(l)((bigz(j)*bigk(i)^(alpha)*l^(1-alpha)...
                    + (1-delta)*bigk(i) - bigk(k))^(-1)*...
                    (1-alpha)*bigz(j)*bigk(i)^(alpha)*l^(-alpha)...
                    + psi*l^(eta-1)); % *U_c(.)F_l(.,.,.) - U_l(.) - Optimal rule for labor
                l(k) = fsolve(foc,1/5,options); %matlab function to solve for the optimal rule for l
                l(k) = real(l(k)); %sometimes fsolve looks for complex number. I want to avoid it!
                l(k) = max(l(k),0.1); %labor should be between zero and 1. I impose nonnegativity
                c(k) = bigz(j)*(bigk(i)^(alpha))*(l(k)^(1-alpha)) + (1-delta)*bigk(i) - bigk(k); 
                %optimal rule for c given l, k, k', and z
                c(k) = max(c(k),0.1); %I want to avoid nonnegativity of consumption
                vnew(k) = log(c(k)) - psi*(l(k)^eta)/eta + beta*v(i,j); %Here, I am defining
                %the dynamic programming problem for given k, k', and z +
                %optimal rule for for labor
            end
            vnew_max(i,j) = max((vnew)); % for every k and z, pick the...
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
