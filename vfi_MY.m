% The planner's problem as VFI, Mendoza & Yue (2012)

betP=bsxfun(@times,bet,P); % Markov probability adjusted discount factor
phiBetP=eye(E)-(1-phi)*betP; % Markov probability and readmission prob adjusted discount factor
rep = repmat(b',B*E,1); % create B*E rows of b'
bp=reshape(rep,B,E*B); % create a B x EB matrix of bond values (tomorrow's bonds)
cL=gdp-L.^omega/omega-bsxfun(@times,xi*log(e),gdp); % "consumption" (2xE, good or bad x shock space)
% (cL is a composite of cons and labor, i.e. it's the input to the utility function)
u_bad=cL(2,:).^(1-sigm)/(1-sigm); % bad utility (1xE)
u_bad=u_bad/phiBetP; % discount it
odf=betP/phiBetP; % overall discount factor of bad scenario
if exist('init_guess','file'), load(f,'D','v'), 
else D=true(B,E); v=zeros(B,E); 
end % load in previous value, else set initial guess to zero for value, ones for decision
errD=1; errV=1; % initialize errors
dvmin=1e-5; % initialize
tic
while errD>0||errV>0
    q=(1-D*P)/r; % set bond price (BxE)
    q(q<0)=0; % make sure it's nonnegative
    D_old=D; errV=1;
    % Below: define utility of good case as good cons + bond*price - b_p[Bx EB]. This is not yet the utility, I'll get to it in a second.
    u_good=repmat(bsxfun(@plus,cL(1,:),bsxfun(@times,q,b)),1,B)-bp; 
    u_good(u_good<0)=0; % make sure it's nonnegative for inputting it to the utility function
    u_good=u_good.^(1-sigm)/(1-sigm); % input it to the utility function. Now it's the utility of the good scenario.
    while errV>dvmin % as long as value has not converged...
        v_bad=u_bad+phi*v(1,:)*odf; % set bad value (1xE)
        v_good=reshape(max(u_good+repmat(v*betP,1,B)),E,B)'; % good value is the max for different b_p choices [BxE]
        D=bsxfun(@gt,v_bad,v_good); % default = 1 if vbad > vgood
        [~,i]=find(D); % find where default occurs (D !=0)
        v_good(D)=v_bad(i); % when you default, replace vg with vb
        errV=max(max(abs(v_good./v-1))); % evaluate the error
        v=v_good; % update value
    end
    errD=sum(D(:)~=D_old(:)); % error = sum the number of times D unequal to D_old
    fprintf('error in D =%g\n',errD)
    if errD==0; dvmin=0; end % if D has converged, set dvmin = 0 to end the big while loop 
end,toc

%Below: I'm defining d_pos to capture both the index of the debt stance
%within state matrix and also the default decision.
% Thus (d_pos = 0 if default,
%       else d_pos = index within debt space b of current debt level)
[~,d_pos]=max(u_good+repmat(v*betP,1,B)); % get the debt position of max values (i.e. index of debt stance where value is maximized)
d_pos=reshape(d_pos,E,B)'; d_pos(D)=0; % set debt position to 0 for default, leave it at whatever other value it was for no default.
save('init_guess','D','v') % save updated v as initial guess
realized_debt=reshape(b(sum(~D)),1,E); % the realized debt choices w/o default
le=log(e); % get log of shock process
E_med = round(E/2); % get median value of shock (for use in graphs and effects at median) 

figure(1),hold on,plot(b,q(:,E_med),'Color','b','linewidth',2), title('Bond price at median TFP'),xlabel('Debt')

figure(2), mesh(Z, b, q), title('Bond price as function of debt and productivity'),xlabel('Debt'), ylabel('log TFP')
