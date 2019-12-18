function [res] = I(T)
tau = T;
b = 0.012;
beta = -0.2;
sigma = 0.14;
k = sqrt( beta.^2 + 2 .* sigma.^2 );
gamma0 = 0.012/0.2;

A = - (2.*b) ./ (sigma .^ 2) .* ...
    log( (2 .* k .* exp( (k - beta) .* tau ./ 2  ) ) ./ ...
        ( (k - beta) .* ( exp( k .* tau ) -1 ) + 2 .* k ) ...
    );
 B = (2 .* ( exp(k.*tau) - 1 ) ) ./ ...
     ( (k - beta) .* ( exp( k.*tau ) - 1 ) + 2 .* k );

 
 
 S = exp( -A -B .* gamma0 );
 res = S * p0(T);
end

