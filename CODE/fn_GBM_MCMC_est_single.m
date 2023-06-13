function out = fn_GBM_MCMC_est_single(in)

% series
Snan   = in.S; 
delta  = in.delta;

% meta
bank = in.bank;
ctry = in.ctry;

% MCMC settings
Niter = in.mcmc.iters;
burn  = in.mcmc.burn;
% setup priors
   mup = in.priors.mu_mean;
  musp = in.priors.mu_stdev;
sigmap = in.priors.sigma_mean;
    v0 = in.priors.v0;

% setup proposals
r          = in.r;
mu_std     = in.mu_std;
sigma_std  = in.sigma_std;
    
% display settings
progr    = in.disp.progr;
plt_post = in.disp.post;
plt_corr = in.disp.corr;
plt_ser  = in.disp.ser;
prt_ests = in.disp.ests;

% deal with NaNs
[S1,mult] = fn_GBM_mult(Snan);

% parameter initialization
   mu = in.init.mu;
sigma = in.init.sigma;

pars = [mu;sigma];

% initialize likelihood
logL = -fn_GBM_loglik([pars(1),log(pars(2))],S1,delta,mult);
% initialize priors
% r-truncated normal prior for mu: N((-4r,r);mup,musp)
%           gamma prior for sigma: G(v0/2,2*sigmap/v0)
logP = log(fn_TNORM_pdf(pars(1),r,mup,musp)) + ... 
       log(gampdf(pars(2),v0/2,2*sigmap/v0)); 

% initialize variables
out_mu     = zeros(Niter,1);
out_sigma  = out_mu;
out_logL   = out_mu;

if plt_ser
    figure
    plot(in.S_old,'*r')
    hold on
    plot(Snan)
    hold off
    title([ctry ': ' bank])
    xlabel('total income')
end

if progr
    h = waitbar(0,'Please wait...');
end
for iter = 1:Niter
    
    if progr
        % progress bar
        waitbar( iter / Niter );   
    end
    
    % store results
    out_mu(iter)    = pars(1);
    out_sigma(iter) = pars(2);
    out_logL(iter)  = - fn_GBM_loglik([pars(1);log(pars(2))],S1,delta,mult);
    
    % mu update
    mu_old = pars(1);
    pars(1) = fn_TNORM_sim(r,pars(1),mu_std,1);
    
    logL1 = -fn_GBM_loglik([pars(1);log(pars(2))],S1,delta,mult);
    logP1 = log(fn_TNORM_pdf(pars(1),r,mup,musp))+...
            log(gampdf(pars(2),v0/2,2*sigmap/v0)); 
    logQ1  = log(fn_TNORM_pdf(mu_old,r,pars(1),mu_std));
    logQ   = log(fn_TNORM_pdf(pars(1),r,mu_old,mu_std)); 
    
    alpha = logL1 - logL + logP1 - logP + logQ1 - logQ;
    if ( alpha >= log( rand ) )
        logL = logL1;
        logP = logP1;
    else
        pars(1) = mu_old;
    end
    
    % sigma update
    M = log(pars(2)/sqrt(1+(sigma_std^2)/(pars(2)^2)));
    V = sqrt(log(1+(sigma_std^2)/(pars(2)^2)));
    sigma_new = M + V*randn(); % log-normal distribution
    logL1 = - fn_GBM_loglik([pars(1),sigma_new],S1,delta,mult);
    logP1 = log(fn_TNORM_pdf(pars(1),r,mup,musp))+...   
            log(gampdf(exp(sigma_new),v0/2,2*sigmap/v0)); 
    logR1 = log(normpdf(log(pars(2)),log(exp(sigma_new)/sqrt(1+(sigma_std^2)/(exp(sigma_new)^2))),...
        sqrt(log(1+(sigma_std^2)/(exp(sigma_new)^2)))));
    logR = log(normpdf(sigma_new,M,V)); 
        
    alpha = logL1 - logL + logP1 - logP + logR1 - logR;
    if ( alpha >= log( rand ) )
        logL = logL1;
        logP = logP1;
        pars(2) = exp(sigma_new);
    end
    
end
if progr
    close(h);
end

if plt_post
    % posterior draws plots
    figure
    subplot(2,1,1)
    plot(out_mu);
    title('Posterior draws of mu');
    ylabel('mu')
    subplot(2,1,2)
    plot(out_sigma);
    title('Posterior draws of sigma');
    ylabel('sigma')
    xlabel('iter')
end

if plt_corr
    % autocorrelation plots
    noLags = 100;
    %
    figure
    subplot(2,1,1)
    autocorr(out_mu,noLags);
    title('autocorrellation of mu')
    ylabel('mu');
    subplot(2,1,2)
    autocorr(out_sigma,noLags);
    title('autocorrellation of sigma')
    ylabel('sigma');
    xlabel('lags');
end

out.mean.mu    = mean(out_mu(burn:end));
out.mean.sigma = mean(out_sigma(burn:end));

out.med.mu     = median(out_mu(burn:end));
out.med.sigma  = median(out_sigma(burn:end));

out.acr.mu     = acrate(out_mu(burn:end));
out.acr.sigma  = acrate(out_sigma(burn:end)); 

out.pror_mu     = fn_TNORM_sim(r,mup,musp,1e5);   
out.pror_sigma  = gamrnd(v0/2,2*sigmap/v0,1e5,1); 

out.post_mu     = out_mu;
out.post_sigma  = out_sigma;

if prt_ests
    % country and bank
    disp(' ')
    disp([ctry ': ' bank])
    disp(' ')
    
    % MCMC estimation results
    disp('Mean point estimates:' );
    disp([{'mu'},{'sigma'}]);
    disp([out.mean.mu,out.mean.sigma]);
    disp('Median point estimates:');
    disp([{'mu'},{'sigma'}]);
    disp([out.med.mu,out.med.sigma]);

    % acceptance rate
    disp('Acceptance rate mu:');
    disp(out.acr.mu);
    disp(' ')
    disp('Acceptance rate sigma:');
    disp(out.acr.sigma);
    
end

end
