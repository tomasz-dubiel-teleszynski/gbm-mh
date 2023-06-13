function out = fn_GBM_MCMC_est_joint(in)

% series
Snan   = in.S; 
delta  = in.delta;

% meta
bank = in.bank;
ctry = in.ctry;

% MCMC settings
Niter  = in.mcmc.iters;
extrad = in.mcmc.iters_extra_adapt;
beta   = in.mcmc.weight_beta; 
disper = in.mcmc.weight_disper;

% setup priors
   mup = in.priors.mu_mean;
  musp = in.priors.mu_stdev;
sigmap = in.priors.sigma_mean;
    v0 = in.priors.v0;
    
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

   pars = [mu;log(sigma)];
   logL = -fn_GBM_loglik(pars,S1,delta,mult);
pars(2) = exp(pars(2));

% initialize priors
% normal prior for mu:    N(mup,0.1)
% gamma prior for sigma: G(v0/2,2*sigmap/v0)

logP = log(normpdf(pars(1),mup,musp))+log(gampdf(pars(2),v0/2,2*sigmap/v0));

% adaptive RW Metropolis settings
d = numel(pars);

% MCMC paramters
burn = floor(Niter/2);
fixpar = 2.38;

% initialize variables
out_mu     = zeros(Niter,1);
out_sigma  = out_mu;
out_logL   = out_mu;

if plt_ser
    figure
    plot(Snan)
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
    
    % update all
    if iter < 2*d + 1 + extrad
         pars1 = fn_mvnrnd([pars(1);log(pars(2))],(disper^2)/d*eye(d),1)';
    else
        covmat = cov([out_mu(1:iter),log(out_sigma(1:iter ))]);
        covmat = (covmat+covmat')/2;
         pars1 = (1-beta)*fn_mvnrnd([pars(1);log(pars(2))],(fixpar^2)/d*covmat,1)'...
            +beta*fn_mvnrnd([pars(1);log(pars(2))],(disper^2)/d*eye(d),1)';
    end
    
    % update likelihood
       logL1 = - fn_GBM_loglik(pars1,S1,delta,mult);
    pars1(2) = exp(pars1(2));
    
    % update prior
    logP1 = log(normpdf(pars(1),mup,musp))+log(gampdf(pars(2),v0/2,2*sigmap/v0));   
 
    % accept-reject
    alpha = logL1-logL+logP-logP1;
    if (alpha >= log(rand))
        pars = pars1;
        logL = logL1;
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
    disp('Acceptance rate:');
    disp(out.acr.mu);
    
end

end
