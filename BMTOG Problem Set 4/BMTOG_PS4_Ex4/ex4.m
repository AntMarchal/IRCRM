clear all;
clc;
% We define the global variables: 
%       b0, b1 and b2 are beta_0 beta_1 and beta_2
%       zcbonds stands for zero coupon bonds.
global maturities b0 b1 b2 a zcbonds
b0 = 0.07;
b1 = -0.02;
b2 = 0.01;
a = 1;
% maturities is a vector containing [0.5, 1, 1.5, ...,9.5,  10]
%   it is used to represent the times of all the possible cash flows
maturities = (1:20)./2;

% zcbonds holds the price of the zero coupon bonds for every maturity
zcbonds = compute_zcoupon_price(maturities);

% we create all the bonds by using the object Bond
%   the input of the constructor is the relevant cash flows of each bond
portfolio = Bond([repmat(10, 1,19) 110])
bond1 = Bond([0 1 repmat(0, 1, 18)])
bond2 = Bond([0 0 0 0 0 1 repmat(0, 1, 14)])
bond3 = Bond([repmat(0,1,9) 1 repmat(0, 1, 10)])

% We compute the weigths of the hedged portfolio by inverting a simple 3x3
% sistem of equations.
A = [[bond1.b0sens bond2.b0sens bond3.b0sens]
     [bond1.b1sens bond2.b1sens bond3.b1sens]
     [bond1.b2sens bond2.b2sens bond3.b2sens]];
 
y = - [[portfolio.b0sens]
       [portfolio.b1sens]
       [portfolio.b2sens]];
   

weights = inv(A)*y

hedged_cash_flows = portfolio.cash_flows ...
+ weights(1)*bond1.cash_flows + weights(2)*bond2.cash_flows ...
+weights(3)*bond3.cash_flows;

% finally we create the hedged bond portfolio. We can check that is has
% zero sensitivity to each beta.
hedged_portfolio = Bond(hedged_cash_flows)