function X = fn_GBM_sim( X0, mu, sigma, delta, T )

N = floor( T / delta );

% set up
X = zeros( N, 1 );
X( 1, 1 ) = log(X0);

% fix values
A = (mu-0.5*sigma^2)*delta;
B = sigma*sqrt(delta);

% simulate signal
for n = 2:N
    X(n,1) = X(n-1,1) + A + B*randn(1,1);
end

X = exp(X);