function [Z,F,B] = tauchen_MY(cover,sige, rho_e, mu, numz, lbb, ubb, numb)

% Tauchen discretization of TFP for Mendoza and Yue
%Preallocating matrices to speed up the code
Z = zeros(numz,1); %discrete spanning space of TFP (I'm calling TFP z)
F = zeros(numz,numz); %Markov chain transition probability matrix

lbz = mu - cover*sige/sqrt((1-rho_e^2)); %lower bound of the z grid
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
B=linspace(lbb,ubb,numb)';