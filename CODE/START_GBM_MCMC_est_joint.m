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
in.init.mu    = init(1);
in.init.sigma = init(2);
% MCMC settings
in.mcmc.iters             = 10000;
in.mcmc.iters_extra_adapt = 0;
in.mcmc.weight_beta       = 0.05;
in.mcmc.weight_disper     = 0.1;
% priors settings
%dev                  =  1.1; % deviation from real values
in.priors.v0         =  5.0;
in.priors.mu_mean    =  0;%dev*mu0;
in.priors.mu_stdev   =  0.5;
in.priors.sigma_mean =  0.2;%dev*sigma0;
% display settings
in.disp.progr = 1;
in.disp.post  = 1;
in.disp.corr  = 1;
in.disp.ests  = 1;
in.disp.ser   = 1;

% run
out = fn_GBM_MCMC_est_joint(in);