function [Call,Put] = Bond_Option_Price(T,S,K,r_0,b,beta,sigma)

%=========================================================================%
%========== Bond Vanilla Option  price for the Vasicek model =============%
%=============== (Note that T,S can be multidimensional) =================%
%=========================================================================%

% Compute a matrix of ZCB prices for the Vasicek Model
% 1st column: P(0,T),  2nd column: P(0,S) 
P = Vasicek_ZCB_Price(0,[T,S],r_0,b,beta,sigma);

% L^2 norm of sigma_T,S
sigma_L_2 = sigma / beta * (exp(beta * (S-T)) - 1)...
         .* sqrt((exp(2 * beta * T) - 1)/(2 * beta));

% Thresholds
d = @(i) (log(P(:,2)./(P(:,1).*K)) - 1/2*(-1)^i * sigma_L_2.^2)./sigma_L_2;

Call = P(:,2) .* normcdf(d(1))  - K .* P(:,1) .* normcdf(d(2))  ;

Put  = K .* P(:,1) .* normcdf(-d(2)) -  P(:,2) .* normcdf(-d(1));

end