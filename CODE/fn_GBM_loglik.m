function logL = fn_GBM_loglik(parameters,S1,delta,mult)

    X = diff(log(S1));
    N = size(X,1);
    
    % unknown paramers
    mu    =     parameters(1);
    sigma = exp(parameters(2));
    
    % define
    sigma2 = sigma^2;
    A = delta*(mu-0.5*sigma2);
    B = delta*sigma2;
    
    logL = 0;
    for t = 1:N
        delta_i = mult(t);
        M = delta_i*A;
        V = delta_i*B;
        logL = logL - 0.5*log(2*pi)-0.5*log(V)...
            - 0.5*((X(t)-M)^2)/V;
    end
    logL = -logL;
    
end
