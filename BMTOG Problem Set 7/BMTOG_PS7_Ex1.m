%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================== Problem Set 7 ============================%
%================================ Exercise 1 =============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format long; warning('off')

%% 0. Setup 

beta = -0.86; b = 0.09 * abs(beta); sigma = 0.0148; r_0 = 0.08;

T_0 = 1/4; T = (T_0:T_0:1)'; n = length(T)-1; c = 4; N = 100;

%% I. Price of the Put Option on the Coupon Bond

% Cash flows of the coupon bond
CF = c +  N * (T == T(end));

% ZCB price
P = @(t,T,r) Vasicek_ZCB_Price(t,T,r,b,beta,sigma);

% Array of coupon dates (excluding T_0)
S = T(2:end);

% Function r -> p(r) (value of the coupon bond at T_0, r(0) = r)
p = @(r) dot(CF(2:end),P(T_0,S,r));

% Strike: T_0-forward price of the coupon bond
K = dot(CF,P(0,T,r_0)) / P(0,T_0,r_0);

% Unique short rate r^* s.t. p(r^*) = K 
r_star = fzero(@(r) p(r)-K,r_0);

fprintf('Value of the short rate r^*: %2.4f\n',r_star)

% Equivalent strikes for a portfolio of ZCB Put Options
K_star = P(T_0,S,r_star);

% Array of ZCB Put Option Prices using the strikes derived above
[~,Put_ZCB] = Bond_Option_Price(repmat(T_0,n,1),S,K_star,r_0,b,beta,sigma);

% Put Option Price of the Coupon Bond
Put_CB = dot(CF(2:end),Put_ZCB);

fprintf('\nPrice of the Put Option on the Coupon Bond: %2.5f\n',Put_CB)
