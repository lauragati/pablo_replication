function filtered_data = hpfilter_lg(data,T,lam)
% data   = is the dataset to be filtered
% T      = length of dataset or simulation length
% lambda = smoothing parameter

% I'm using the Matlab helpfile of hpfilter to recreate this command.
hp=repmat([lam -4*lam 1+6*lam -4*lam lam],T,1);
hp([T+1 2*T-1 3*T+2 4*T])=-2*lam; 
hp([2*T+1 3*T])=1+lam; 
hp([2*T+2 3*T-1])=1+5*lam;
HP=spdiags(hp,-2:2,T,T);
filtered_data = data-HP\data;