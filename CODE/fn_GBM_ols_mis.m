function out = fn_GBM_ols_mis(in)

Snan   = in.S; 
delta  = in.delta;

[S1,mult] = fn_GBM_mult(Snan);

dlogS = diff(log(S1))./sqrt(mult*delta);
Y = dlogS;
N = length(dlogS);
X = sqrt(mult*delta);

b = (X'*X)\X'*Y;
e = Y-X*b;

df = numel(b)+1;

vOLS = 1/(N-df)*(e'*e);
vMLE = (N-df)/N*vOLS;

sigma = sqrt(vMLE);
mu = b+0.5*vMLE;

out.mu = mu;
out.sigma = sigma;

end
