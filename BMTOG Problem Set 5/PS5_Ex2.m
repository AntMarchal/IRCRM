%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================== Problem Set 5 ============================%
%================================ Exercise 2 =============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format short; warning('off')

%% 0. Setup 

beta = -0.86; b = 0.09 * abs(beta); sigma = 0.0148; r_0 = 0.095;

T = 2; N_T = T * 252; D_t = T / N_T; N_sim = 10000;

grid = linspace(0,T,N_T+1);

%% I. Monte Carlo Exact Scheme

% Store the terms exp(beta*t_i) in a vector 
Y = exp(beta * grid');

% Exact mean and variance of the process
Mu = r_0 * Y + b/beta * (Y-1); Var = sigma^2*(Y.^2-1)/(2*beta);

% Resulting trajectories for the short rate 
r = Mu + sqrt(Var).* [ones(1,N_sim);randn(N_T,N_sim)];

%% II. Evolution of the Bank Account and obtained price

Bank_Account = [ones(1,N_sim);exp(cumsum(r(1:end-1,:))/252)];

Terminal_Value = mean(Bank_Account(end,:));

MC_Price = 1 / Terminal_Value;

fprintf('\nMonte Carlo Price of the ZCB Bond: %2.5f\n',MC_Price)

%% III. Plots

figure 

subplot(2,1,1); plot(grid,r(:,1:200),'Linewidth',1.5); 
xlabel('t');    title('Short Rate')

subplot(2,1,2); plot(grid,Bank_Account(:,1:200),'Linewidth',1.5);
xlabel('t');    title('Bank Account')

%% Comparison with the closed form solution and exact moments

y = Y(end); B = (y - 1)/beta;

A = sigma^2*(4*y-y^2-2*log(y)-3)/(4*beta^2) + b*(y-1-log(y))/(beta^2);

fprintf('\nExact Price of the ZCB Bond: %2.5f\n',exp(-(A + B * r_0)))

% Mean
fprintf('\n====================================================\n')
fprintf('\nMonte Carlo mean of r at time T: %2.5f\n',mean(r(end,:)))

mean_exact = r_0 * exp(beta * T) + b/beta * (exp(beta*T)-1);

fprintf('\nExact mean of r at time T: %2.5f\n',mean_exact)

% Variance
fprintf('\n====================================================\n')
fprintf('\nMonte Carlo variance of r at time T: %e\n',var(r(end,:)))

var_exact = sigma^2/(2*beta) * (exp(2*beta*T)-1);

fprintf('\nExact variance of r at time T: %e\n',var_exact)
