function x = fn_TNORM_sim(r,mu,sigma,N)

mlt = 4;
% truncate from (-mlt*r,r)
x = mu + sigma * norminv( rand(N,1)*(normcdf((r-mu)/sigma,0,1)-normcdf((-mlt*r-mu)/sigma,0,1))+normcdf((-mlt*r-mu)/sigma,0,1) ,0,1);

end

