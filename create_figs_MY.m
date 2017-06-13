% Recreate Table III and Figure VI using the solution of VFI and simulating the model

%% Simulation
T=500; % total number of periods in simulation including burn-in
Ts=400; % number of periods in simulations minus burn-in
smooth=1600; %the lambda smoothing parameter for HP filter
n_sim=2000; % number of simulations
num_var=8; % number of variables
sim_means=nan(n_sim,num_var); % means of simulated variables
sim_var=nan(n_sim,num_var);  % variances of simulated variables
sim_corr=nan(n_sim,7); % correlations of simulated variables
historical=nan(n_sim,3); % called 'historical' because it refers to the 'historical default-output comovement' from MY Table III.
Y_def=[]; % all variables around (equal periods before and after) default; for all simulations
debtmax=12; %debt upper bound levels for "before-and-after default" window
debt_window=-debtmax:debtmax; % debt levels for "before-and-after default" window
F=cumsum(P); %take a cumulative Markov probability in order to determine Markov chain of TFP shocks
for sim=1:n_sim
    Y=nan(T,num_var); % variables vector
    def_hist=false(T,1); %default history, initialized all as false (no default)
    reentry_hist=def_hist; %reentry history, initialzed all as false (no reentry) - will be set in a second.
    H=nan(T,3); % default, reentry & state history gathered in one matrix
    d=false; i=1; j=E_med; % initial states (index of default, debt and TFP shock respectively)
    tfp_index=rand(T,1); % draw TFP indexes, they will be used to compute the chain of TFP draws at end of each loop
    reentry=rand(T,1)<phi; % draw reentry for each simulation
    for t=1:T
        running_index=d_pos(i,j); % index of default for this period
        if running_index==0||d % if index is zero, i.e. there is default or if there was default yesterday, there's default today
            if ~d % if yesterday there was no default, change state d to default since today we're defaulting.
                d=true; 
                def_hist(t)=d; % and add it to default history
            end 
            % also, set all variables to their default values...
            R=nan; % sovereign bond spread is nan in default b/c no interest rate exists on 0 bonds
            GDP=gdp(2,j); % real GDP
            TB=xi*log(e(j))*GDP; % trade balance
            MS=ms(2,j);% imported intermediate inputs
            IM=M(2,j); % total intermediate goods
            L_=L(2,j); % labor (total)
            running_index=i; %... and update the running index.
        else %If there's no default today, set variables to their good values
            R=1/q(running_index,j)^4-r^4; % sovereign bond spread (annual because MY use annual)
            GDP=gdp(1,j); 
            TB=b(i)-b(running_index)*q(running_index,j)+xi*log(e(j))*GDP; 
            MS=ms(1,j); 
            IM=M(1,j); 
            L_=L(1,j); 
        end
        CONS=GDP-TB; % consumption
        Y(t,:)=[GDP CONS R TB/GDP L_ IM MS b(i)/GDP]; % gather all variables in one vector
        H(t,:)=[d i j]; % record the current state in the history vector
        % For next period:
        if d && reentry(t) % If we can reenter markets tomorrow, then...
            d=false; % ...set the default indicator for false (so that tomorrow it's possible not to default)
            reentry_hist(t)=true; % ...add the reentry to credit markets to reentry history
            running_index=1;  %...and change running index. 
        end 
        j=find(tfp_index(t)<F(:,j),1); % update TFP index to get draw for tomorrow using Markov chain
        i=running_index; %update the debt position of tomorrow according to how we updated the running index (which depended on whether we default or not, and whether we reenter or not)
    end
    Y(:,[1:2 6:7]) =log(Y(:,[1:2 6:7])); %take logs for variables for which it makes sense
    Y(:,[1 2 4 6 7])=hpfilter_lg(Y(:,[1 2 4 6 7]),T,smooth); % invoke HP filter (MY use HP to detrend so I'm following them)
    %Below: construct the matrix of variable values before and after default.
    default_index=bsxfun(@plus,find(def_hist(T-Ts+1:T))+T-Ts,debt_window)'; %get the time periods around default (as indexes of a window), and add them to the debt level 
    Y_def_sim=nan(size(default_index,1)*size(default_index,2),num_var+2); % initialize vector to be size of default occurances x number of variables plus 2 b/c we'll add the history of default and reentry too
    di=default_index>=1&default_index<=T; %get all the indexes that are between 1 and T
    default_index=default_index(di);
    Y_def_sim(di,:)=[Y(default_index,:) def_hist(default_index) reentry_hist(default_index)]; % vector of variables around times of default (in a window before and after); for this simulation only
    Y_def=[Y_def; Y_def_sim]; % in each loop, I add the variables of this simulation, so this matrix expands in each loop.
    %after the loop, the size of Y_def is (number of default occurences) x (number of variables)
    sim_means(sim,:)=nanmean(Y(T-Ts+1:T,:)); % calculate means of variables in this simulation 
    sim_var(sim,:)=nanvar(Y(T-Ts+1:T,:)); % calculate variances of variables in this simulation 
    CC=nancov(Y(T-Ts+1:T,1:6))./sqrt(sim_var(sim,1:6)'*sim_var(sim,1:6));
    sim_corr(sim,:)=[CC(1,3:6) CC(3,4:6)]; % calculate correlation coefficients of variables in this simulation 
    CC=corrcoef([Y(T-Ts+1:T,1) def_hist(T-Ts+1:T)]);
    autocorr_GDP=corrcoef(Y(T-Ts+1:T-1,1),Y(T-Ts+2:T,1)); % autocorr GDP
    historical(sim,:)=[CC(2) mean(def_hist(H(:,2)>1)) autocorr_GDP(2)];
end
Y_def=permute(reshape(Y_def,2*debtmax+1,size(Y_def,1)/(2*debtmax+1),num_var+2),[1 3 2]); % reshape this matrix so that it now is (number of debt levels) x (number of variables) x (number of defaults per each debt level) --> so that I can take means across all simulations
Y_def_means=nanmean(Y_def,3); %take means along 3rd dimension, i.e.  we get means for all variables across all simulations around times of default
%% Table III
Y=mean(sim_means); % take means of simulated means (for all simulations)
sim_var=sqrt(mean(sim_var)); % calculate standard deviations of simulated variables at the mean (for all simulations)
sim_corr=mean(sim_corr); % calculate means of standard deviations of simulated variables (for all simulations)
historical=nanmean(historical); 
fprintf('Average debt/GDP ratio %.2f%%\n',Y(8)*100) 
fprintf('Average bond spreads %.2f%%\n',Y(3)*100) 
fprintf('Std. dev. of bond spreads %.2f%%\n',sim_var(3)*100) 
fprintf('Consumption std.dev./GDP std.dev. %.2f\n',sim_var(2)/sim_var(1))
fprintf('Correlations with GDP\n')
fprintf('   bond spreads %.2f\n',sim_corr(1))
fprintf('   trade balance %.2f\n',sim_corr(2))
fprintf('   labor %.2f\n',sim_corr(3))
fprintf('   intermediate goods %.2f\n',sim_corr(4)) 
fprintf('Correlations with bond spreads\n') 
fprintf('   trade balance %.2f\n',sim_corr(5)) 
fprintf('   labor %.2f\n',sim_corr(6)) 
fprintf('   intermediate goods %.2f\n',sim_corr(7)) 
fprintf('Historical default-output comovement\n')
fprintf('   correlation between default and GDP %.2f\n',historical(1)) 
fprintf('   fraction of defaults with GDP below trend %.2f%%\n',mean(Y_def(debtmax+1,1,:)<0)*100) 
fprintf('   fraction of defaults with large recessions %.2f%%\n',mean(Y_def(debtmax+1,1,:)<-2*sim_var(1))*100) %where GDP is less than 2 variances below trend

%% Figure VI
figure(3),titles={'GDP','Consumption','Interest rate','Trade balance/GDP','Labor','Intermediate goods','Imported intermediate goods','Debt/GDP'};
for i=1:num_var
    subplot(2,4,i),hold on
    plot(debt_window,Y_def_means(:,i),'r-','LineWidth',2,'Color','k')
    set(gca,'XTick',-debtmax:4:debtmax)
    set(gca,'XTickLabel',-debtmax/4:debtmax/4)
    title(titles{i}),xlabel('year')
    axis([-debtmax debtmax -inf inf]),grid on
end