function [y,P,s]=rouwenhorst(rho,sigma,n)
% Rowenhurst method to approximate univariate process by Markov chain
%      y(t) = rho y(t-1)+ e(t)
% where y(t) is AR(1) process
%       e(t) is iid random variable with mean 0 and std.dev. sigma
% INPUTS: rho - autocorrelation coefficient
%         sigma - std.dev. of e(t)
%         n - number of states in Markov chain
% Reference:
% "Finite State Markov-Chain Approximations to Highly Persistent Processes."
% by Karen A. Kopecky and Richard M. H. Suen.
% Review of Economic Dynamics 13 (2010), pp. 701-714.
% URL: http://www.karenkopecky.net/RouwenhorstPaperFinal.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: Iskander Karibzhanov
%          Department of Economics
%          University of Minnesota
%          karib003@umn.edu

ybar=sqrt((n-1)/(1-rho^2))*sigma;
y=linspace(-ybar,ybar,n);
p=(1+rho)/2; q=p;
P=rh(n);
s=zeros(n,1);
for j=1:n
    s(j)=nchoosek(n-1,j-1);
end
s=s/2^(n-1);

    function P=rh(h)
        if h==2
            P=[p 1-p; 1-q q];
        else
            P1=rh(h-1);
            z=zeros(1,h);
            z1=zeros(h-1,1);
            P=[p*P1 z1; z]+[z1 (1-p)*P1; z]+...
              [z; (1-q)*P1 z1]+[z; z1 q*P1];
            P(2:h-1,:)=P(2:h-1,:)/2;
        end
    end
end