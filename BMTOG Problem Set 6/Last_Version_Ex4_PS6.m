%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================== Problem Set 6 ============================%
%================================ Exercise 4 =============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

%% Setup
close all; clear; clc; format long; warning('off') 

%% Parameters
beta = -0.5; b = 0.5*0.02; sigma = 0.01; r_0 = 0.01;
dt= 0.25; 
T = 30;
T_0 = 0.25;
n=T/dt;

% array of time
t   = T_0:dt:T;

%% Model
% Computation of B, A and P
for i = 1:length(t)
    B = 1/beta * ( exp(beta*t(i) - 1));
    A = sigma^2 /(4*beta^3) * (4*exp(beta*t(i))...
         - exp(2*beta*t(i))-2*beta*t(i)-3)...
         + b * (exp(beta*t(i))-1-beta*t(i))/beta^2;
    P(i) = exp(-A-B*r_0);
end

% Rswap
Rswap = (P(1) - P(end)) / (dt*sum(P(2:end)));

% Price of each caplet
for i = 2:length(t)

    caplet2(i-1) = caplet(dt,beta,sigma,P(i),P(i-1),Rswap,t(i-1));
end

% Price of Cap
cap = sum(caplet2)

% Archive Code
function capletPrice = caplet(dt,beta,sigma,P,P_,R,T)

d1 = (log((1 + dt * R)*P/P_)...
   + sigma^2/(4*beta^3)*(exp(beta*dt)-1)^2*(exp(2*beta*T)-1))...
   / sigma/beta*(exp(beta*dt)-1)*sqrt((exp(2*beta*T)-1)/(2*beta));

d2 = (log((1 + dt * R)*P/P_)...
   - sigma^2/(4*beta^3)*(exp(beta*dt)-1)^2*(exp(2*beta*T)-1))...
   / sigma/beta*(exp(beta*dt)-1)*sqrt((exp(2*beta*T)-1)/(2*beta));

capletPrice =  P_*(1-normcdf(d2)) - (1 + dt * R)* P *(1-normcdf(d1));

end
