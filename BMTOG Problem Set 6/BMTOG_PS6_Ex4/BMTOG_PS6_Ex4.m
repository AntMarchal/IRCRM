%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================== Problem Set 6 ============================%
%================================ Exercise 4 =============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format long; warning('off')


% %% 0. Setup 
% 
% beta = -0.5; b = 0.5*0.02; sigma = 0.01; r_0 = 0.01;
% dt= 0.25; T=30; T0=0;
% % T = 30*252;
% % T_0 = 252/4;
% mat=[0:dt:T]';
% %%
% % 
% 
% n= T/dt-1;
% 
% for i = 2:n
%     B(i) = 1/beta * ( exp(beta*mat(i)) - 1);
%     A(i) = 1/(4*beta^3) * (sigma^2 *(4*exp(beta*mat(i))-exp(2*beta*mat(i))-2*beta*mat(i)-3));
%     P(i) = exp(-A(i)-B(i)*r_0);
% end
% 
% R_swap = (P(n) - P(1)) / (dt*sum(P));
% for i=2:n
%     capletPrice(i)=caplet(dt,beta,sigma,P(i),P(i-1),R_swap,mat(i));
% end
% 
% cap=sum(capletPrice)
% 
% %% Functions
% 
% 
% function capletPrice = caplet(dt,beta,sigma,P,P_,r,T)
% 
% d1=(log(P/(r*P_))+1/(2*beta)*sigma*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/(2*beta)))...
%     /( sigma/beta*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/2*beta));
% d2 = (log(P/(r*P_))-1/(2*beta)*sigma*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/(2*beta)))...
%     /( sigma/beta*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/2*beta));
% 
% capletPrice =  (1/dt-r)*P_*(1-normcdf(d2)) - P/dt*(1-normcdf(d1));
% end


%%
%% Setup
close all; clear; clc; format long; warning('off') 
%% Parameters
beta = -0.5; b = 0.5*0.02; sigma = 0.01; r_0 = 0.01;
dt= 0.25; 
T = 30;
T_0 = 0.25;
% array of time
mat = zeros(120,1);
mat(1) = T_0;
for i = 2:120
    mat(i) = mat(i-1) + T_0;
end
%% Model
% Computation of B, A and P
for i = 1:120
    B(i) = 1/beta * ( exp(beta*mat(i)) - 1);
    A(i) = 1/(4*beta^3) * (sigma^2 *(4*exp(beta*mat(i))-exp(2*beta*mat(i))-2*beta*mat(i)-3));
    P(i) = exp(-A(i)-B(i)*r_0);
end
% Rswap
Rswap = (P(1) - P(end)) / (dt*sum(P));

% L
for i = 1:120
    L(i) = (1/mat(i)) * ( (1/P(i)) - 1);
end

% Price of each caplet
for i = 1:120
    if i>1
    caplet(i) = (1/P(i)) * (mat(i) - mat(i-1)) * max(0, (L(i) - Rswap));
    else
    caplet(i) =(1/P(i)) * (mat(i)) * max(0, (L(i) - Rswap))
    end
end

% Price of Cap
cap = sum(caplet)
