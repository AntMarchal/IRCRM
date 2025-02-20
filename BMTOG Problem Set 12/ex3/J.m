function [res] = J(s)
b = 0.012;
beta = -0.2;
sigma = 0.14;
k = sqrt( beta.^2 + 2 .* sigma.^2 );
tau = s;
gamma0 = 0.012/0.2; % not given in the problem set

A = - (2.*b) ./ (sigma .^ 2) .* ...
    log( (2 .* k .* exp( (k - beta) .* tau ./ 2  ) ) ./ ...
        ( (k - beta) .* ( exp( k .* tau ) -1 ) + 2 .* k ) ...
    );
 B = (2 .* ( exp(k.*tau) - 1 ) ) ./ ...
     ( (k - beta) .* ( exp( k.*tau ) - 1 ) + 2 .* k );

 
 dA = -(2 .* b) ./ (sigma.^2) .* ...
     ( ...
            (k - beta) ./ (2)...
             -  ...
             ( (k -beta) .* ( exp(k.*tau) - 1 ) + 2 .* k ) .^ (-1) .* ( k - beta ) .* k .* exp( k .* tau ) ...
     );
 
 dB = B .* (...
         ( 2 .* ( exp( k.*tau ) - 1 ) ).^(-1)  .* 2 .* k .* exp( k .* tau ) ...
         - ... 
          ( ( k - beta ) .* ( exp( k .* tau ) - 1 ) + 2 .* k ).^(-1) .* ...
          ( k - beta ) .* k .* exp( k .* tau )...
     );
 
 S = exp( -A -B .* gamma0 );
 res =  S .* (dA + dB .* gamma0) .* p0(s);
end

