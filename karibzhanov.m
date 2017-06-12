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
fun=@(L,ps,ej) [ps*gamma*A*L(2)^(gamma-1)*(gamma*alpha_m/alpha_L*L(1)/L(2)-1)^(1/mu-1)*(lambda/(1-lambda))^(1/mu)-sum(L)^(omega-1)
     ej*alpha_L*L(1)^(alpha_L-1+alpha_m/mu)*(A*L(2)^(gamma-1/mu)*(alpha_m/alpha_L*gamma*lambda)^(1/mu))^alpha_m-sum(L)^(omega-1)];
ps=(theta*[Rstar^(nu/(nu-1));0]+1-theta).^(1-1/nu);
L=nan(2,E,2); % regime(g|b) x schock(e) x sector(f|m)
o=optimset('Display','off');
for i=1:2
    for j=1:E
        L(i,j,:)=fsolve(@(L) fun(L,ps(i),e(j)),[.07 .05],o);
    end
end
Lf=L(:,:,1); Lm=L(:,:,2); L=sum(L,3);
w=L.^(omega-1);
y=w.*Lf/alpha_L; % final production
md=A*Lm.^gamma;
pm=w.*Lm./md/gamma;
ms=md.*bsxfun(@rdivide,pm,ps/(1/lambda-1)).^(1/(1-mu));
if mu==0; M=md.^lambda.*ms.^(1-lambda);
else M=(lambda*md.^mu+(1-lambda)*ms.^mu).^(1/mu); end
gdp=y-bsxfun(@times,ps,ms); %applies the element-by-element binary operation specified by the function handle FUNC to arrays A and B
fprintf('Model %d\n',model)
disp(['GDP autocorrelation (no defaults): ',num2str((pd.*(gdp(1,:)-gdp(1,:)*pd'))*P'*gdp(1,:)'./var(gdp(1,:),pd))])
disp(['GDP std.dev. (no defaults): ',num2str(sqrt(var(gdp(1,:),pd))*100),'%'])
figure(1),hold on
plot(log(e(ei)),100*(1-gdp(2,ei)./gdp(1,ei)),[color '.-']),title('Output drop in default as a function of TFP Shock'),xlabel('log of TFP shock'),ylabel('GDP drop, %')
line(log(e([1 E])),100*(1-gdp(2,[1 E])./gdp(1,[1 E])),'LineStyle','--','Color','k')
disp('Effects of default on factor allocations')
fprintf('M\t%.2f\n',(M(2,E0)/M(1,E0)-1)*100)
fprintf('m*\t%.2f\n',(ms(2,E0)/ms(1,E0)-1)*100)
fprintf('md\t%.2f\n',(md(2,E0)/md(1,E0)-1)*100)
fprintf('L\t%.2f\n',(L(2,E0)/L(1,E0)-1)*100)
fprintf('Lf\t%.2f\n',(Lf(2,E0)/Lf(1,E0)-1)*100)
fprintf('Lm\t%.2f\n',(Lm(2,E0)/Lm(1,E0)-1)*100)
%% solve dynamic problem
A=501;%501
amax=.05; a=linspace(0,amax,A)'; % discretize bond space
g=exp(Gamma+zg(2,:));
bP=bsxfun(@times,beta*g.^(1-sigma),P); iP=eye(E)-(1-phi)*bP;
aa=reshape(repmat(a',A*E,1),A,E*A);
cL=gdp-L.^omega/omega-bsxfun(@times,xi*log(e),gdp);
ub=cL(2,:).^(1-sigma)/(1-sigma); ub=ub/iP; bPi=bP/iP;
if exist(f,'file'), load(f,'D','v'), else D=true(A,E); v=zeros(A,E); end
dD=1; dv=1; dvmin=1e-5; tic
figure(2),hold on,title('Default set'),xlabel('Debt'),ylabel('log of TFP schock')
while dD>0||dv>0
    q=(1-D*P)/Rstar; q(q<0)=0; D_old=D; dv=1;
    ug=repmat(bsxfun(@plus,cL(1,:),bsxfun(@times,q,a*g)),1,A)-aa;
    ug(ug<0)=0; ug=ug.^(1-sigma)/(1-sigma);
    while dv>dvmin
        vb=ub+phi*v(1,:)*bPi;
        vg=reshape(max(ug+repmat(v*bP,1,A)),E,A)';
        D=bsxfun(@gt,vb,vg); [~,i]=find(D); vg(D)=vb(i);
        dv=max(max(abs(vg./v-1))); v=vg;
    end
%     plot(reshape(a(sum(~D))+a(2)/2,Z,G)',reshape(log(e),Z,G)',color),drawnow
    dD=sum(D(:)~=D_old(:));
    fprintf('dD=%g\n',dD)
    if dD==0; dvmin=0; end
end,toc
[~,ap]=max(ug+repmat(v*bP,1,A)); ap=reshape(ap,E,A)'; ap(D)=0; save(f,'D','v') % save updated v as initial guess
ad=reshape(a(sum(~D))+a(2)/2,Z,G); le=reshape(log(e),Z,G);
figure(2),hold on,plot(ad,le,color)
title('Default set'),xlabel('Debt'),ylabel('log of TFP shock')
figure(3),hold on,plot(a,q(:,[1 E]),'Color',color), title('Bond price'),xlabel('Debt')
%% Simulation
T=500; Ts=400; % number of periods in simulations
s=1600; S=repmat([s -4*s 1+6*s -4*s s],T,1);
S([T+1 2*T-1 3*T+2 4*T])=-2*s; S([2*T+1 3*T])=1+s; S([2*T+2 3*T-1])=1+5*s;
HP=spdiags(S,-2:2,T,T); % HP filter
nSim=2000; Sn=8; Sm=nan(nSim,Sn); Sv=nan(nSim,Sn); Sr=nan(nSim,7); STAT=nan(nSim,3); Sd=[];
dtmax=12; dt=-dtmax:dtmax; F=cumsum(P);
for sim=1:nSim
    S=nan(T,Sn); Dh=false(T,1); Rh=Dh; H=nan(T,3); % default, rehab & state histroy
    h=false; i=1; j=E0; % initial states
    s=rand(T,1); re=rand(T,1)<phi; % TFP shocks, rehab draws
    for t=1:T
        i1=ap(i,j); % index of default for this period
        if i1==0||h % if index is zero, i.e. there is default or if there was default yesterday
            if ~h; h=true; Dh(t)=h; end % if yesterday there was no default, change that.
            R=nan;
            GDP=gdp(2,j); 
            TB=xi*log(e(j))*GDP;
%             MD=md(2,j);
            MS=ms(2,j);
            IM=M(2,j);
            L_=L(2,j);
            i1=i;
        else
            R=1/q(i1,j)^4-Rstar^4; % annual sovereign bond spread
            GDP=gdp(1,j); % real GDP
            TB=a(i)-g(j)*a(i1)*q(i1,j)+xi*log(e(j))*GDP; % trade balance
%             MD=md(1,j); % domestic intermediate inputs
            MS=ms(1,j); % imported intermediate inputs
            IM=M(1,j); % total intermediate goods
            L_=L(1,j); % labor
        end
        CON=GDP-TB; % consumption
        S(t,:)=[GDP CON R TB/GDP L_ IM MS a(i)/GDP];
        H(t,:)=[h i j];
        % next period:
        if h && re(t); h=false; Rh(t)=true; i1=1; end
        j=find(s(t)<F(:,j),1);
        i=i1;
    end
%     if any(Dh),figure(4),plot(ad,le,'-',a(H(:,2)),log(e(H(:,3))),'.',a(H(find(H(:,1),1)-(-15:2),2)),log(e(H(find(H(:,1),1)-(-15:2),3))),'r.-'),pause, end;
    gh=[0;cumsum(Gamma+zg(2,H(1:T-1,3)))'];
    S(:,1:2)=log(S(:,1:2))+[gh gh];
    S(:,6:7)=log(S(:,6:7));
    S(:,[1 2 4 6 7])=S(:,[1 2 4 6 7])-HP\S(:,[1 2 4 6 7]); % H-P detrended
%     S(:,6:7)=detrend(log(S(:,6:7))); % log-linearly detrended
    d=bsxfun(@plus,find(Dh(T-Ts+1:T))+T-Ts,dt)';
    sd=nan(size(d,1)*size(d,2),Sn+2);
    di=d>=1&d<=T; d=d(di);
    sd(di,:)=[S(d,:) Dh(d) Rh(d)];
    Sd=[Sd; sd]; %#ok<AGROW>
    Sm(sim,:)=nanmean(S(T-Ts+1:T,:));
    Sv(sim,:)=nanvar(S(T-Ts+1:T,:));
    CC=nancov(S(T-Ts+1:T,1:6))./sqrt(Sv(sim,1:6)'*Sv(sim,1:6));
    Sr(sim,:)=[CC(1,3:6) CC(3,4:6)];
    CC=corrcoef([S(T-Ts+1:T,1) Dh(T-Ts+1:T)]);
    AC=corrcoef(S(T-Ts+1:T-1,1),S(T-Ts+2:T,1));
    STAT(sim,:)=[CC(2) mean(Dh(H(:,2)>1)) AC(2)];
end
Sd=permute(reshape(Sd,2*dtmax+1,size(Sd,1)/(2*dtmax+1),Sn+2),[1 3 2]);
Sd(:,5,:)=bsxfun(@rdivide,Sd(:,5,:),Sd(1,5,:));
Sdm=nanmean(Sd,3);
Sdv=nanstd(Sd,0,3);
Sdm_ce=nanmean(Sd(:,1:Sn,all(Sd(dtmax+1:2*dtmax+1,Sn+2,:)==0)),3); % continued exclusion
Sdm_re=nanmean(Sd(:,1:Sn,Sd(dtmax+1,Sn+2,:)==1),3); % immediate reentry
%% Statistical Moments
S=mean(Sm); Sv=sqrt(mean(Sv)); Sr=mean(Sr); STAT=nanmean(STAT);
fprintf('Average debt/GDP ratio %.2f%%\n',S(8)*100) %
fprintf('Average bond spreads %.2f%%\n',S(3)*100) %
fprintf('GDP autocorrelation %.2f\n',STAT(3))
fprintf('Std. dev. of GDP %.2f%%\n',Sv(1)*100)
fprintf('Std. dev. of trade balance %.2f%%\n',Sv(4)*100)
fprintf('Std. dev. of bond spreads %.2f%%\n',Sv(3)*100) %
fprintf('Consumption std.dev./GDP std.dev. %.2f\n',Sv(2)/Sv(1))%
fprintf('Correlations with GDP\n')
fprintf('   bond spreads %.2f\n',Sr(1))%
fprintf('   trade balance %.2f\n',Sr(2))%
fprintf('   labor %.2f\n',Sr(3))%
fprintf('   intermediate goods %.2f\n',Sr(4)) %
fprintf('Correlations with bond spreads\n') 
fprintf('   trade balance %.2f\n',Sr(5)) %
fprintf('   labor %.2f\n',Sr(6)) %
fprintf('   intermediate goods %.2f\n',Sr(7)) %
fprintf('Historical default-output co-movements\n')
fprintf('   correlation between default and GDP %.2f\n',STAT(1)) %
fprintf('   fraction of defaults with GDP below trend %.2f%%\n',mean(Sd(dtmax+1,1,:)<0)*100) %
fprintf('   fraction of defaults with large recessions %.2f%%\n',mean(Sd(dtmax+1,1,:)<-2*Sv(1))*100) %
fprintf('Quaterly default frequency %.2f%%\n',STAT(2)*100)
fprintf('Output drop in default %.2f%%\n',(Sdm(dtmax+1,1)-Sdm(dtmax,1))*100)
if ~exist('stat.txt','file'), f = fopen('stat.txt','w'); fprintf(f,'Debt\tSpread\tGDPac\tGDPstd\tTBstd\tRstd\tCstd\tR&GDP\tTB&GDP\tL&GDP\tIM&GDP\tTB&R\tL&R\tIM&R\tD&GDP\tDlow\tDrec\tDfreq\tGDPdrop\n');
else f = fopen('stat.txt','a'); end
fprintf(f,'%.2f\t',[S(8)*100 S(3)*100 STAT(3) Sv(1)*100 Sv(4)*100 Sv(3)*100 Sv(2)/Sv(1) Sr STAT(1) mean(Sd(dtmax+1,1,:)<0)*100 mean(Sd(dtmax+1,1,:)<-2*Sv(1))*100 STAT(2)*100 (Sdm(dtmax+1,1)-Sdm(dtmax,1))*100]);
fprintf(f,'\n'); fclose(f);
%% macro dynamics around default events
figure(4),titles={'GDP','Consumption','Interest rate','Trade balance/GDP','Labor','Intermediate goods','Imported intermediate goods','Debt/GDP'};
for i=1:Sn
    subplot(2,4,i),hold on
    plot(dt,Sdm(:,i),'r-','LineWidth',2,'Color',color)
%    plot(dt,[Sdm(:,i) Sdm_re(:,i) Sdm_ce(:,i)],'-','LineWidth',2),hold on
%    plot(dt,[Sdm(:,i)+Sdv(:,i),Sdm(:,i)-Sdv(:,i)],'k--')
    set(gca,'XTick',-dtmax:4:dtmax)
    set(gca,'XTickLabel',-dtmax/4:dtmax/4)
    title(titles{i}),xlabel('year')
    axis([-dtmax dtmax -inf inf]),grid on
end
%h={'Model average','Immediate reentry','Continued exclusion'};
%legend(h,'Orientation','horizontal','Position',[0.4 0.5 0.24 0.025]);