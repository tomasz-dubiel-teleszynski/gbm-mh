function [S1,mult,idx] = fn_GBM_mult(S)
% get varying steps
idx = ~isnan(S);
S1 = S(idx);
mult = diff(find(idx));
idx = find(idx);
end