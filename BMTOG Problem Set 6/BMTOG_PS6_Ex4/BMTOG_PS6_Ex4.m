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
mat=[0.5:dt:T];
%% Model
% Computation of B, A and P
for i = 1:119
    B(i) = 1/beta * ( exp(beta*dt - 1));
    A(i) = 1/(4*beta^3) * (sigma^2 *(4*exp(beta*dt)-exp(2*beta*dt)-2*beta*dt-3))...
        +b*( exp(beta*dt)-1-beta*dt)/beta^2;
    P(i) = exp(-A(i)-B(i)*r_0);
end

% Rswap
Rswap = (P(1) - P(end)) / (dt*sum(P(2:end)));
% L
for i = 1:119
    L(i) = (1/mat(i)) * ( (1/P(i)) - 1);
end

% Price of each caplet
for i = 2:119
    caplet2(i) = caplet(dt,beta,sigma,P(i),P(i-1),L(i),mat(i));
end

% Price of Cap
cap = sum(caplet2)

% Archive Code
function capletPrice = caplet(dt,beta,sigma,P,P_,r,T)

d1=(log(P/(r*P_))+1/(2*beta)*sigma*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/(2*beta)))...
    /( sigma/beta*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/2*beta));
d2 =(log(P/(r*P_))-1/(2*beta)*sigma*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/(2*beta)))...
    /( sigma/beta*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/2*beta));

capletPrice =  (1/dt-r)*P_*(1-normcdf(d2)) - P/dt*(1-normcdf(d1));
end
