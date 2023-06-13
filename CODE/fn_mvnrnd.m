function x = fn_mvnrnd(mu,sig2,n)

[k, l] = size(mu);
if k>l
    mu = mu';
end

sig2 = (sig2+sig2')/2+1e-6*eye(numel(mu));

m = numel(mu);
L = chol(sig2);
x = repmat(mu,n,1) + randn(n,m)*L;

end

