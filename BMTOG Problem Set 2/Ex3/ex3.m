clc
%% Part a)
% We define the following global variables:
%   maturities: a vector of the maturity dates [1, 2, ... 10]
%   yields: a vector of all yields corresponding to each maturity date:
%       [0.06, 0.058, ... , 0.051]
%   zcbonds: a vector containing the price of the zero coupon bonds 
%   hedged_cash_flows: a vector containing the cash flows of the 
%       hedged portfolio

global maturities yields zcbonds hedged_cash_flows
maturities = transpose(1:10);
yields = [6, 5.8, 5.62, 5.46, 5.33, 5.25, 5.2, 5.16, 5.125, 5.1].'/100;
zcbonds = exp( - maturities .* yields );

% We created an object called Bond that represents the bonds
%   We initialize three bonds. The initialization computes automatically
%   the price, duration and convexity of the bonds and prints them.
   
portfolio = Bond([6, 8, 106, 7, 0, 102, 3, 3, 3, 110].')
bond1 = Bond([4, 4, 4, 4, 4, 4, 104].')
bond2 = Bond([10,10,10,10,10,10,10,110].')


%% Part b)


% Implement a convexity hedge:
%       We compute the weights of the convexity hedge by
%       solving a 2 equation problem in matrix form.

A=[(bond1.duration*bond1.price) (bond1.convexity*bond1.price);
   (bond2.duration*bond2.price) (bond2.convexity*bond2.price)]
weights = ( [  (-portfolio.price*portfolio.duration) (-portfolio.price*portfolio.convexity)]) * inv(A)

%% Part c)
% We Plot the evolution of the hedged portfolio with shifts of the yield curve

% We define the interval for s:
X= (-0.2):0.001:0.2;

% We compute the cash flows of the hedged portfolio as a linear combination
% of the three bond's cash flows
hedged_cash_flows = portfolio.cash_flows + [ bond1.cash_flows ; [0 0 0].'] * weights(1) + [ bond2.cash_flows; [0 0].' ]* weights(2);

% Finally we create the hedged_cash_flows bond. We can observe that it's 
% duration and convexity are infinitesmal
Bond(hedged_cash_flows)

% We define the difference in the value of the portfolio
% hedged_of_s(s) computes the value of the hedged portfolio with a shift 
% of s in the yield curve.
Y= hedged_of_s(X)- hedged_of_s(0);
plot(X,Y)
xlabel('shift')
ylabel('Value change')
