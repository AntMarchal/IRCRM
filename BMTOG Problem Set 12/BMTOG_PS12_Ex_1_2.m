%=========================================================================%
%================= Interest Rate and Credit Risk Models ==================% 
%============================ Problem Set 12 =============================%
%=========================================================================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format long; warning('off')

%% Exercise 1
fprintf("============= Exercise 1 ==============\n")

% Parameters (d:delta, D:Delta_t,l:lambda)
r = .05; d = 0.6; D = 1/2; T = 10; l = .02;

%============================= Part (b) ==================================%

% Fair CDS spread
x_star = @(l,d,r)  d/D * l /(r+l) * (exp((r+l) * D) - 1);

x_star_1 = x_star(l,d,r);

fprintf("\nFair CDS spread: %2.2f bps\n",1e4 * x_star_1)

%============================= Part (c) ==================================%

% Notional Amount ($)
N = 1e8;

V_prem = @(x) x * D * (1- exp(-(r+ l) * T))/(exp((r+l) * D) - 1);

Pi_Y = N * (V_prem(x_star_1) - V_prem(1e-2));

fprintf("\nAmount to pay upfront: $%2.1f\n", Pi_Y)

%% Exercise 2

% Parameters (d:delta)
r = .01; d = 0.4; T = 1;

%============================= Part (a) ==================================%

% Fair CDS spread (half year and one year respectively)
x_0 = .02; x_1 = .04;

% Compute a = lambda / 2
a = fzero(@(l) x_star(l,d,r) - x_0,0) /2;

fprintf("\nParameter a: %2.4f\n",a)

%============================= Part (b) ==================================%

% Premium leg as a function of b
V_prem = @(b) x_1/2 * exp(-(r + 2*a)/2) * (1 + exp(-(r + 2*a + b)/2));

% Default leg as a function of b
V_def  = @(b) d * (2*a/(r + 2*a) * (1 - exp(-(r + 2*a)/2))...
       + (2*a+b)*exp(-(r + 2*a)/2)/(r+2*a+b)*(1 - exp(-(r + 2*a + b)/2)));
   
b = fzero(@(b) V_prem(b) - V_def(b),0);

fprintf("\nParameter b: %2.4f\n",b)
  



