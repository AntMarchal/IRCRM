function [res] = c(T, sigma)
global L r;
VbyB = 1./L;
dp = ( log(VbyB) + ( r + sigma.^2./2 ) .* T ) ./ (sigma .* T.^0.5 );
dm = ( log(VbyB) + ( r - sigma.^2./2 ) .* T ) ./ (sigma .* T.^0.5 );
res = -1./T .* log( 1 + exp(r.* T) .* VbyB .* normcdf(-dp) - normcdf(-dm) );
end


