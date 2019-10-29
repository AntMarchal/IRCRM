%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================== Problem Set 6 ============================%
%================================ Exercise 4 =============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format long; warning('off')


%% 0. Setup 

beta = -0.5; b = 0.5*0.02; sigma = 0.01; r_0 = 0.01;
dt= 0.25; T=30; T0=0;
T = 30*252;
T_0 = 252/4;
t(1) = T_0;
n = T/252 * 4;
for i = 2:n
    t(i) = t(i-1) + T_0;
end

for i = 1:n
    B(i) = 1/beta * ( exp(beta*t(i)) - 1);
    A(i) = 1/(4*beta^3) * ( sigma^2 * (4*exp(beta*t(i)) - exp(2*beta*t(i)) - 2 * beta * t(i) - 3) );
    P(i) = exp(-A(i) - B(i)*r_0);
end
R_swap = (P(n) - P(1)) / (T_0*sum(P));


for i = 1:n
    L(i) = (1/t(i)) * ( (1/P(i)) - 1);
end

for i = 1:n
    if i>1
    f(i) = 1/dt * (P(i)/P(i-1) - 1);
    else
    f(i) = 1/dt * (1/P(i) - 1);
    end
    
end

for i = 1:n
    if i>1
    capletPrice(i) = caplet(R_swap, f(i), sigma, L(i), t(i), t(i-1));
    else
    capletPrice(i) = caplet(R_swap, f(i), sigma, L(i), t(i), 0);
    end
end

cap=sum(capletPrice)


%% Setup

%%
close all; clear; clc; format long; warning('off')

k = 0.5;
theta = 0.02;
r0 = 0.01;
sigma = 0.01;
days = 252;
years = 30;
T = years*days;
T0 = days/4;
beta = -k;
b = theta * k;
n = T/days * 4;

t = zeros(n,1);

t(1) = T0;
for i = 2:n
    t(i) = t(i-1) + T0;
end

P = zeros(n,1);
A = zeros(n,1);
B = zeros(n,1);

for i = 1:n
    B(i) = -1/k * ( exp(-k*t(i)) - 1);
    A(i) = 1/(4*beta^3) * ( sigma^2 * (4*exp(beta*t(i)) - exp(2*beta*t(i)) - 2 * beta * t(i) - 3) );
    P(i) = exp(-A(i) - B(i)*r0);
end

R_swap = (P(1) - P(end)) / (T0*sum(P));

delta = t(2) - t(1);

L = zeros(n,1);
forward = zeros(n,1);
forward(1) = 1/delta * (1/P(1) - 1);

for m = 1:n
    L(m) = (1/t(m)) * ( (1/P(m)) - 1);
end

for j = 2:n
    forward(j) = 1/delta * (P(i)/P(i-1) - 1);
end

cplt = zeros(n,1);

cplt(1) = cpltprice(R_swap, forward(1), sigma, L(1), t(1), 0);

for i = 2:n
    cplt(i) = cpltprice(R_swap, forward(i), sigma, L(i), t(i), t(i-1));
end

cap = sum(cplt)

function cplprice = cpltprice(strike, fwdrate, vol, spot_rate, t, t_)
d1 = real( ( log(fwdrate/strike) + vol^2 * t_/2 ) / (vol * sqrt(t_) ));
d2 = real(d1 - vol * sqrt(t_));
cplprice = (fwdrate * normcdf(d1) - strike * normcdf(d2)) / (1+spot_rate)^t;
end

function capletPrice = caplet(k, f, sigma, r, t, t_)
d1 = real( ( log(f/k) + sigma^2 * t_/2 ) / (sigma * sqrt(t_) ));
d2 = real(d1 - sigma * sqrt(t_));
capletPrice = (f * normcdf(d1) - k * normcdf(d2)) / (1+r)^t;
end
