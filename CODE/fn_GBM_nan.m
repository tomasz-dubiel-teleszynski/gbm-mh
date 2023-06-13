function S = fn_GBM_nan(S)
% 25% of time series is NaN
N = length(S);
cap = ceil(N/4); 
idx = randi(N,N,1);
S(idx<cap) = NaN;

end