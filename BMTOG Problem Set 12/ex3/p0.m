function [res,A,B] = p0(T)

global r;
tau = T;
b = 0.077;
beta = -0.96;
sigma = 0.014;

A = sigma.^2 .* ( 4 .* exp( beta .* tau ) - exp( 2 .* beta .* tau ) - 2 .* beta .* tau - 3 ) ./ ...
    (4 .* beta.^3) + ...
    b .* ( exp( beta .* tau ) - 1 - beta .* tau ) ./ ...
    ( beta.^2 );
 B = 1 ./ beta .*...
     ( exp(beta .* tau) - 1 );
 S = exp( -A -B .* r );
 res = S;
end

