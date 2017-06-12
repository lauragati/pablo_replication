clear all; model=2; color='r'; % 1: only transitory shocks, 2: only growth shocks, 3: both
f=['model' num2str(model) '.mat'];% rng(1)
% if exist(f,'file'), delete(f), end
% close all; clc; %if exist('stat.txt','file'), delete('stat.txt'), end
alpha_m=.43; % share of intermediate goods in gross output
gamma=.7; % labor income share in value added and in intermediate goods production
alpha_L=(1-alpha_m)*gamma;
alpha_k=1-alpha_m-alpha_L;
sigma=2; % coef. of relative risk aversion
Rstar=1.01; % risk-free interest rate
omega=1.455; % curvature of labor supply, 1/(omega-1)=Frisch wage elasticity
phi=.083; % reentry probability
lambda=.62; % armington weight of domestic inputs
mu=.65; % armington curvature parameter
nu=.59; % Dixit-Stiglitz curvature parameter
A=.31; % intermediate goods TFP coefficient
beta=.88; % subjective discount factor
theta=.7; % upper bound of imported inputs with working capital
xi=0;%-.67; % TFP semi-elasticity of exogenous capital flows
% autocorrelation, std.dev. and mean of transitory & growth shocks:
Gamma=.02; %sigma_e=.05; wg=1; E=101; % weight on g=[0,1]
Z=1; rho_z=.95; sigma_z=0; mu_z=-sigma_z^2/(1-rho_z^2)/2; 
G=101; rho_g=0;   sigma_g=.017; mu_g=-sigma_g^2/(1-rho_g^2)/2;
% Z=ceil(E^(1-wg)); rho_z=.95; sigma_z=(1-wg)*sigma_e*sqrt(1-rho_z^2); mu_z=-sigma_z^2/(1-rho_z^2)/2; 
% G=ceil(E^wg); rho_g=0; sigma_g=wg*sigma_e*sqrt(1-rho_g^2); mu_g=-sigma_g^2/(1-rho_g^2)/2;
% if model==1, G=1; elseif model==2, Z=1; end
E=Z*G; E0=round(E/2); Z0=round(Z/2); G0=round(G/2);
% [zg,P,pd]=rouwenhorst(rho_z,sigma_z,Z); zg=[zg+mu_z;zeros(1,Z)]; P=P'; pd=pd';
% zg=[linspace(-1,1,Z)*(sigma_z*sqrt(3/(1-rho_z^2)*(Z-1)/(Z+1)))+mu_z; zeros(1,Z)];
% P=repmat((1-rho_z)/Z,Z,Z); P(1:Z+1:Z*Z)=1-(Z-1)/Z*(1-rho_z); pd=repmat(1/Z,1,Z);
[zg,P,pd]=tauchen([Z G],[mu_z mu_g],diag([rho_z rho_g]),[sigma_z sigma_g],[20 20]);
e=exp(sum(zg)); [~,ei]=sort(e); % ei becomes the index of e (~ is being ignored by def.)
%% solve static problem
fun=@(L,ps,ej) [alpha_m*ej*(lambda*(A*L(2)^gamma)^mu + (1-lambda)*L(3)^mu)^((alpha_m-mu)/mu)*L(1)^alpha_L*(1-lambda)*L(3)^(mu-1) - ps
    alpha_m*ej*(lambda*(A*L(2)^gamma)^mu + (1-lambda)*L(3)^mu)^((alpha_m-mu)/mu)*L(1)^alpha_L * lambda*(A*L(2)^gamma)^(mu-1) - (L(1)+L(2))^(omega-1)/(gamma*A*L(2)^(gamma-1))
    alpha_L*ej*(lambda*(A*L(2)^gamma)^mu + (1-lambda)*L(3)^mu)^(alpha_m/mu)*L(1)^(alpha_L -1) - (L(1)+L(2))^(omega-1)];
ps=(theta*[Rstar^(nu/(nu-1));0]+1-theta).^(1-1/nu);
L=nan(2,E,3); % regime(g|b) x schock(e) x sector(f|m) and ms
o=optimset('Display','off');
for i=1:2
    for j=1:E
        L(i,j,:)=fsolve(@(L) fun(L,ps(i),e(j)),[.07 .05 .005],o);
    end
end
Lf=L(:,:,1); Lm=L(:,:,2); L=sum(L,3);

% it works and yields same solution as that of Karibzhanov! So adopt this
% version for vfi_v4