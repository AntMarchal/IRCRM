%=========================================================================%
%================= Interest Rate and Credit Risk Models ==================% 
%============================ Problem Set 11 =============================%
%============================== Exercise 2 ===============================%
%=========================================================================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format short; 

%% 0. Setup

% State Space
S = ["A","B","C","D"]; 

% Ratings of firms 1,2,3
R = ["A","B","C"]; 

% Initial asset value
V_0 = 100;

% Transition Matrix
P = [.95 .05 0 0;.05 .8 .1 .05;0 .2 .5 .3;0 0 0 1];

% Drift and volatility of the firms
mu_V  = [.2 .15 .1]; sig_V = [.3 .25 .2];

%% I. Thresholds of Transition

% Summary table
Thresholds = table("d_tilde_" + string(4:-1:0)','VariableNames',"Threshold");
                       
% Function matching rating to corresponding vector index
ind = @(J) find(ismember(S,J));

% Standardized thresholds for firms 1,2,3
d = norminv(cumsum(P(ind(R),:),2,'reverse'))';

% Transform the thresholds back into the original scale
d_tilde = [V_0 * exp((mu_V - 1/2*sig_V.^2) + sig_V.*d); zeros(1,3)];

% Display the results
Thresholds = [Thresholds,array2table(d_tilde,'VariableNames',R)]
