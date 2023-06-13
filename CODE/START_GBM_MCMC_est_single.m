clear all
close all
clc

% set real values
     S0 = 1;
    mu0 = 0.01;
 sigma0 = 0.05;
      T = 20;
  delta = 0.25;
% simulate GBM process
   S = fn_GBM_sim(S0,mu0,sigma0,delta,T);
Snan = fn_GBM_nan(S);
  S1 = fn_GBM_mult(Snan);
init = fn_GBM_ols(S1,delta);

in = struct();
% meta
in.bank = 'bank';
in.ctry = 'ctry';
% series
in.S     = Snan;
in.delta = delta;
% initialization
in.init.mu    = 0;   %init(1);
in.init.sigma = 0.2; %init(2);
% MCMC settings
in.mcmc.iters = 10000;
in.mcmc.burn  = floor(0.1*in.mcmc.iters);  
% proposals settings
in.r          = 0.06; % truncation
in.mu_std     = 0.4;
in.sigma_std  = 0.125;
% priors settings
in.priors.mu_mean    =  0;
in.priors.mu_stdev   =  0.1;
in.priors.sigma_mean =  0.2;
in.priors.v0         =  5.0;
% display settings
in.disp.progr = 1;
in.disp.post  = 1;
in.disp.corr  = 1;
in.disp.ests  = 1;
in.disp.ser   = 1;

% outliers via prior predictive
in = fn_GBM_prior_pred(in);
% run
out = fn_GBM_MCMC_est_single(in);

% plot distributions
[a_mu_pror,b_mu_pror]=ksdensity(out.pror_mu);
[a_mu_post,b_mu_post]=ksdensity(out.post_mu);
[a_sigma_pror,b_sigma_pror]=ksdensity(out.pror_sigma);
[a_sigma_post,b_sigma_post]=ksdensity(out.post_sigma);

figure
subplot(2,1,1)
plot(b_mu_pror,a_mu_pror,'b')
hold on
plot(b_mu_post,a_mu_post,'r')
hold off
title('mu')
xlabel('prior (blue), posterior (red)')
subplot(2,1,2)
plot(b_sigma_pror,a_sigma_pror,'b')
hold on
plot(b_sigma_post,a_sigma_post,'r')
hold off
title('sigma')
xlabel('prior (blue), posterior (red)')
