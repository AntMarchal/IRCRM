%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================ Problem Set 10 ==============================%
%=========================================================================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format long; warning('off')

%% Part b)

%================================ Setup ==================================%
Ratings = {'A','B','C'}; T_max = 20;

% Table containing the default probabilities
Default_Probs = table('Size',[T_max,4],...
                      'VariableTypes',repmat("double",1,4),...
                      'VariableNames',['T',Ratings]);

Default_Probs.T = (1:T_max)';
                  
% Transition matrix over one time step
P = 1e-2 * [95,5,0,0;5,80,10,5;0,20,50,30;0,0,0,100]; 

% Transition matrix over several time steps (gradually updated)
P_T = P;

%=============================== Loop ====================================%
for T = 1:T_max

    % Compute the probabilities of default
    Default_Probs{T,2:end} = (eye(3,4) * P_T * ((1:4)'==4))';
    
    if T < T_max   
        % Transition matrix over T+1 time steps
        P_T = P_T * P;
    end
end

Default_Probs

%============================= Plot ======================================%
figure 
plot(Default_Probs.T, Default_Probs{:,2:end},'Linewidth',1.5); xlim([1,T_max])
xlabel("Time horizons"); ylabel("Default probabilities")
legend(Ratings,'Location','Best'); title("Default probabilities")

%% Part c)

% Market Data
CR = readtable('CR.csv','ReadRowNames' ,true,'ReadVariableNames' ,true);

% Remove empty rows
CR = CR(1:20,:);

% State space 
S = [Ratings,'D']; 

% Estimated transition matrix
P_hat = zeros(4,4);

for j = 1:4

    % Boolean matrix of companies with rating j, for t = 0,...,18
    B_j = ismember(CR{1:end-1,:},S(j));
    
    % Number of companies that are rated j at time t, t = 0,...,18
    N_j = sum(B_j,2);
    
    % Estimate the j-th row of P hat
    for k = 1:4
        
        % Boolean matrix of companies with rating k, for t = 1,...,19
        B_k = ismember(CR{2:end,:},S(k));
        
        % Number of companies rated j at t and k at t + 1, t = 0,...,18
        N_jk = sum(B_j .* B_k,2);
        
        % Maximum Likelihood estimator of p_jk
        P_hat(j,k) = sum(N_jk) / sum(N_j);
             
    end
end

P_hat