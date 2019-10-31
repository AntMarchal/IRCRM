function [P] = Vasicek_ZCB_Price(t,T,r,b,beta,sigma)

%=========================================================================%
%================= ZCB price for the Vasicek model =======================%
%============= (Note that T can be an array of maturities) ===============%
%=========================================================================%

% Introduce a temporary variable to simplify the expressions for A and B
y = exp(beta * (T-t));

B = (y - 1)/beta;

A = sigma^2 * (4*y - y.^2 - 2*log(y) - 3)/(4*beta^3) + (B-(T-t))*b/beta;

% ZCB Price(s)
P = exp(- A - B * r);


