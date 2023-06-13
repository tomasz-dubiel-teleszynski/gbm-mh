function [mu,sigma] = fn_crude_mu_sigma(S1,Snan,delta)

N = find(~isnan(Snan),1,'last')-find(~isnan(Snan),1,'first');
mu = (S1(end)/S1(1))^(1/(N*delta))-1;
sigma = std(diff(log(S1))/delta);

end

