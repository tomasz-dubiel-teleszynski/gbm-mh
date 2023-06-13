function F = fn_TNORM_pdf(x,r,mu,sigma)

mlt = 4;
% truncate from (-mlt*r,r)
F = ( normpdf((x-mu)/sigma,0,1) ) / ( normcdf( (r-mu)/sigma ,0,1) - normcdf((-mlt*r-mu)/sigma,0,1) ) / sigma;

end

