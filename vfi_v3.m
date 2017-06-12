%% A 3rd version of VFI for Mendoza and Yue, inspired by Karibzhanov's code

clear all; 
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
rho_z=.95; sigma_z=0.017; % sigma_z = 0.017
mu_z=0;%-sigma_z^2/(1-rho_z^2)/2; 

numz = 101; % (25) number of TFP realizations for Tauchen's method of TFP discretization
numb = 11; % number of bonds in bond space. Not given by MY, so I set to 11, so the 6th one is 0.
cover = 3; % coverage parameter for Tauchen. This value MY don't tell, so I set to default value of 3
lbb = -10; % lower bound for bond space
ubb = 10; % uppe bound for bond space

[Z,P,~] = tauchen_MY(cover,sigma_z, rho_z, mu_z, numz, lbb, ubb, numb); 
e=exp(Z'); % shock process
[~,ei]=sort(e); % ei becomes the index of e (~ is being ignored by def.)
E = size(e,2); % =numz 
%% Factor market equilibrium
% Simplify 8-dimensional eq system to a 2-dimensional one in Lf and Lm
fun=@(L,ps,ej) [ps*gamma*A*L(2)^(gamma-1)*(gamma*alpha_m/alpha_L*L(1)/L(2)-1)^(1/mu-1)*(lambda/(1-lambda))^(1/mu)-sum(L)^(omega-1)
     ej*alpha_L*L(1)^(alpha_L-1+alpha_m/mu)*(A*L(2)^(gamma-1/mu)*(alpha_m/alpha_L*gamma*lambda)^(1/mu))^alpha_m-sum(L)^(omega-1)];
ps=(theta*[Rstar^(nu/(nu-1));0]+1-theta).^(1-1/nu);
L=nan(2,E,2); % regime(good|bad) x schock(e) x sector(final|intermediate)
o=optimset('Display','off');
for i=1:2
    for j=1:E
        L(i,j,:)=fsolve(@(L) fun(L,ps(i),e(j)),[.07 .05],o);
    end
end
% then solve the other variables residually:
Lf=L(:,:,1);
Lm=L(:,:,2);
L=sum(L,3);
w=L.^(omega-1);
y=w.*Lf/alpha_L; % final production
md=A*Lm.^gamma;
pm=w.*Lm./md/gamma;
ms=md.*bsxfun(@rdivide,pm,ps/(1/lambda-1)).^(1/(1-mu));
M=(lambda*md.^mu+(1-lambda)*ms.^mu).^(1/mu);
gdp=y-bsxfun(@times,ps,ms); %applies the element-by-element binary operation specified by the function handle 

E0 = round(E/2); % look at the effect in the middle of shock space
disp('Effects of default on factor allocations in percent:')
fprintf('M\t%.2f\n',(M(2,E0)/M(1,E0)-1)*100)
fprintf('m*\t%.2f\n',(ms(2,E0)/ms(1,E0)-1)*100)
fprintf('md\t%.2f\n',(md(2,E0)/md(1,E0)-1)*100)
fprintf('L\t%.2f\n',(L(2,E0)/L(1,E0)-1)*100)
fprintf('Lf\t%.2f\n',(Lf(2,E0)/Lf(1,E0)-1)*100)
fprintf('Lm\t%.2f\n',(Lm(2,E0)/Lm(1,E0)-1)*100)

%% The planner's problem as VFI
B=501; % 501
bmax=.05; b=linspace(0,bmax,B)'; % discretize bond space
bP=bsxfun(@times,beta,P); % Markov probability adjusted discount factor
iP=eye(E)-(1-phi)*bP; % Markov probability and readmission prob adjusted discount factor
rep = repmat(b',B*E,1); % create B*E rows of b'
bp=reshape(rep,B,E*B); % create a B x EB matrix of bond values (tomorrow's bonds)
cL=gdp-L.^omega/omega-bsxfun(@times,xi*log(e),gdp); % "consumption" (2xE, good or bad x shock space)
% (cL is a composite of cons and labor, i.e. it's the input to the utility function)
ub=cL(2,:).^(1-sigma)/(1-sigma); % bad utility (1xE)
ub=ub/iP; % discount it
bPi=bP/iP; % overall discount factor of bad scenario
if exist('init_guess','file'), load(f,'D','v'), else D=true(B,E); v=zeros(B,E); end 
% load in previous value, else set initial guess to zero for value, ones for decision
dD=1; dv=1; % initialize errors
dvmin=1e-5; % initialize
tic
while dD>0||dv>0
    q=(1-D*P)/Rstar; % set bond price (BxE)
    q(q<0)=0; % make sure it's nonnegative
    D_old=D; dv=1;
    % define utility of good case as good cons + bond*price - b_p [B x EB]
    ug=repmat(bsxfun(@plus,cL(1,:),bsxfun(@times,q,b)),1,B)-bp; % deleted a *g from a here
    ug(ug<0)=0; ug=ug.^(1-sigma)/(1-sigma);
    while dv>dvmin % as long as value has not converged...
        vb=ub+phi*v(1,:)*bPi; % set bad value (1xE)
        vg=reshape(max(ug+repmat(v*bP,1,B)),E,B)'; % good value is the max for different b_p choices [BxE]
        D=bsxfun(@gt,vb,vg); % default = 1 if vb > vg
        [~,i]=find(D); % find where default occurs (D !=0)
        vg(D)=vb(i); % when you default, replace vg with vb
        dv=max(max(abs(vg./v-1))); % evaluate the error
        v=vg; % update value
    end
    dD=sum(D(:)~=D_old(:)); % error = sum the number of times D unequal to D_old
    fprintf('dD=%g\n',dD)
    if dD==0; dvmin=0; end % if D has converged, set dvmin = 0 to end the big while loop 
end,toc
[~,ap]=max(ug+repmat(v*bP,1,B)); % get indices of max values
ap=reshape(ap,E,B)'; ap(D)=0; % set indices to 0 for default, 1 for no default.
save('init_guess','D','v') % save updated v as initial guess
ad=reshape(b(sum(~D))+b(2)/2,1,E); % the realized debt choices w/o default
%le=reshape(log(e),Z,G);
le=log(e);
figure(1), plot(ad,le,'r','linewidth',2)
title('Default set'),xlabel('Debt'),ylabel('log TFP')

figure(2),hold on,plot(b,q(:,[E0]),'Color','b','linewidth',2), title('Bond price at median TFP'),xlabel('Debt')
% takes 3 min
figure(3), mesh(Z, b, q), title('Bond price as function of debt and productivity'),xlabel('Debt'), ylabel('log TFP')