function X = fn_GBM_sim_pred( X0, mu, sigma, delta)

N = size(mu,1);

% set up
X = log(X0);

% fix values
A = (mu-0.5*sigma.^2)*delta;
B = sigma*sqrt(delta);

% simulate signal
X = X + A + B.*randn(N,1);

X = exp(X);