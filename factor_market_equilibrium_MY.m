% Factor market equilibrium of Mendoza & Yue (2012)

% Simplify 8-dimensional eq system to a 3-dimensional one in Lf, Lm and ms
fun=@(L,ps,ej) [alpha_m*ej*(lam*(A*L(2)^gam)^mu + (1-lam)*L(3)^mu)^((alpha_m-mu)/mu)*L(1)^alpha_L*(1-lam)*L(3)^(mu-1) - ps
    alpha_m*ej*(lam*(A*L(2)^gam)^mu + (1-lam)*L(3)^mu)^((alpha_m-mu)/mu)*L(1)^alpha_L * lam*(A*L(2)^gam)^(mu-1) - (L(1)+L(2))^(omega-1)/(gam*A*L(2)^(gam-1))
    alpha_L*ej*(lam*(A*L(2)^gam)^mu + (1-lam)*L(3)^mu)^(alpha_m/mu)*L(1)^(alpha_L -1) - (L(1)+L(2))^(omega-1)];
ps=(thet*[r^(nu/(nu-1));0]+1-thet).^(1-1/nu);
% Below: initialize solution:
L=nan(2,E,3); % Size is for regime(good|bad) x shock(e) x {sector(final|intermediate) plus ms}
o=optimset('Display','off');
for i=1:2
    for j=1:E
        L(i,j,:)=fsolve(@(L) fun(L,ps(i),e(j)),[.07 .05 .005],o); % solve for Lf, Lm and ms for each TFP realization, and for the good and bad scenario.
    end
end
Lf=L(:,:,1); Lm=L(:,:,2); ms=L(:,:,3); % Separate out solution.
L= Lf + Lm; % define total labor.
w=L.^(omega-1); % wage.
y=w.*Lf/alpha_L; % final production
md=A*Lm.^gam; % domestic intermediate
pm=w.*Lm./md/gam; % price of imported intermediate
M=(lam*md.^mu+(1-lam)*ms.^mu).^(1/mu); %Armington aggregate of the two intermediates
gdp=y-bsxfun(@times,ps,ms); % GDP.
