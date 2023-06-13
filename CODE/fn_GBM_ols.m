function start_values = fn_GBM_ols(S,delta)

dlogS = diff(log(S));
N = length(dlogS);

X = ones(N,1);
Y = dlogS;
b = (X'*X)\X'*Y;
e = Y-X*b;
df = numel(b)+1;
vOLS = 1/(N-df)*(e'*e);
vMLE = (N-df)/N*vOLS;

sigma = sqrt(vMLE/delta);
mu = b/delta+0.5*sigma^2;

start_values = [mu;sigma];

end

