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
T = 30; dt= 0.25;

%%
n= T/dt;
r=nan(n+1,1);
r(1,:)=r_0;
for i=2:n+1
    dr = (b+beta*r(i-1))*dt+sigma*(sqrt(dt)*randn(1));
    r(i)=r(i-1)+dr;
end

mat = [0:dt:T]';

B= 1/beta * ( exp(beta*mat) - 1);
A= 1/(4*beta^3) * (sigma^2 *(4*exp(beta*mat)-exp(2*beta*mat)-2*beta*mat-3));

P = exp(-A-B*r_0);
R_swap = (P(1) - P(end)) / (T*sum(P));
for i=2:n+1
    capletPrice(i)=caplet(dt,beta,sigma,P(i),P(i-1),r(i),mat(i));
end

cap=n*sum(capletPrice)


%% Functions


function capletPrice = caplet(dt,beta,sigma,P,P_,r,T)

d1=real((log(P/(r*P_))+1/(2*beta)*sigma*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/(2*beta)))...
    /( sigma/beta*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/2*beta)));
d2 = real((log(P/(r*P_))-1/(2*beta)*sigma*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/(2*beta)))...
    /( sigma/beta*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/2*beta)));

capletPrice =  (1/dt-r)*P_*(1-normcdf(d2)) - P/dt*(1-normcdf(d1));
end
