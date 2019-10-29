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


n= T/dt-1;
t=[0:dt:30]';

%%
% 
for i = 1:n
    B(i) = 1/beta * ( exp(beta*t(i)) - 1);
    A(i) = 1/(4*beta^3) * (sigma^2 *(4*exp(beta*t(i))-exp(2*beta*t(i))-2*beta*t(i)-3));
    P(i) = exp(-A(i)-B(i)*r_0);
end

R_swap = (P(n) - P(1)) / (dt*sum(P));

for i = 1:n
    L(i) = (1/t(i)) * ( (1/P(i)) - 1);
end

for i = 1:n
    if i>1
    capletPrice(i) = (1/P(i)) * (t(i) - t(i-1)) * max(0, (L(i) - R_swap));
    else
    capletPrice(i) = (1/P(i)) * (t(i)) * max(0, (L(i) - R_swap));
    end
end

cap=sum(capletPrice)

%% Functions

function capletPrice = caplet(k, f, sigma, r, t, t_)
d1 = real( ( log(f/k) + sigma^2 * t_/2 ) / (sigma * sqrt(t_) ));
d2 = real(d1 - sigma * sqrt(t_));
capletPrice = (f * normcdf(d1) - k * normcdf(d2)) / (1+r)^t;
end

% function capletPrice = caplet(dt,beta,sigma,P,P_,r,T)
% 
% d1=(log(P/(r*P_))+1/(2*beta)*sigma*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/(2*beta)))...
%     /( sigma/beta*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/2*beta));
% d2 = (log(P/(r*P_))-1/(2*beta)*sigma*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/(2*beta)))...
%     /( sigma/beta*(exp(beta*dt)-1)*sqrt(exp(2*beta*(T-1)-1)/2*beta));
% 
% capletPrice =  (1/dt-r)*P_*(1-normcdf(d2)) - P/dt*(1-normcdf(d1));
% end





% %%
% %% Setup
% close all; clear; clc; format long; warning('off') 
% %% Parameters
% beta = -0.5; b = 0.5*0.02; sigma = 0.01; r_0 = 0.01;
% dt= 0.25; 
% T = 30;
% T_0 = 0.25;
% % array of time
% mat = zeros(120,1);
% mat(1) = T_0;
% for i = 2:120
%     mat(i) = mat(i-1) + T_0;
% end
% %% Model
% % Computation of B, A and P
% for i = 1:120
%     B(i) = 1/beta * ( exp(beta*mat(i)) - 1);
%     A(i) = 1/(4*beta^3) * (sigma^2 *(4*exp(beta*mat(i))-exp(2*beta*mat(i))-2*beta*mat(i)-3));
%     P(i) = exp(-A(i)-B(i)*r_0);
% end
% 
% % Rswap
% Rswap = (P(1) - P(end)) / (dt*sum(P));
% % L
% for i = 1:120
%     L(i) = (1/mat(i)) * ( (1/P(i)) - 1);
% end
% 
% % Price of each caplet
% for i = 1:120
%     if i>1
%     caplet(i) = (1/P(i)) * (mat(i) - mat(i-1)) * max(0, (L(i) - Rswap));
%     else
%     caplet(i) =(1/P(i)) * (mat(i)) * max(0, (L(i) - Rswap))
%     end
% end
% 
% % Price of Cap
% cap = sum(caplet)
