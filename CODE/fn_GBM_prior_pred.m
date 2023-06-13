function in = fn_GBM_prior_pred(in)


[S1,mult,idx] = fn_GBM_mult(in.S);

N = length(S1)-1;

pred_mu    = fn_TNORM_sim(in.r,in.priors.mu_mean,in.priors.mu_stdev,1e5);
pred_sigma = gamrnd(in.priors.v0/2,2*in.priors.sigma_mean/in.priors.v0,1e5,1);

bands95 = zeros(N,2);
outliers = zeros(N,1);
stepback = 0;
for n = 1:N
    if stepback == 0
        X = fn_GBM_sim_pred(S1(n), pred_mu, pred_sigma, mult(n)*in.delta);
        bands95(n,:) = [quantile(X,0.025), quantile(X,0.975)];
        if ~(S1(n+1)>bands95(n,1) && S1(n+1)<bands95(n,2))
           outliers(n) = 1;
           stepback = stepback + 1;
        end
    else
        X = fn_GBM_sim_pred(S1(n-stepback), pred_mu, pred_sigma, sum(mult(n-stepback:n))*in.delta);
        bands95(n,:) = [quantile(X,0.025), quantile(X,0.975)];
        if ~(S1(n+1)>bands95(n,1) && S1(n+1)<bands95(n,2))
           outliers(n) = 1;
           stepback = stepback + 1;
        else
           stepback = 0;
        end
    end
end

outliers = logical([0;outliers]);
in.S_old = in.S;
in.S(idx(outliers))=NaN;

