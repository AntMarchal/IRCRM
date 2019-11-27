%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================ Problem Set 9 ==============================%
%=========================================================================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format short; warning('off')

%% Market Data

% Implied volatilities
Imp_Vol = 1e-2*[29.3,29.3,29.3,20.8,20.8,18.3,18.3,17.8,17.8,...
                16.3,16.3,16.7,16.7,16.1,16.1,15.7,15.7,15.7,15.7]';
              
% Forward LIBOR curve         
L_0 = 1e-2*[4.228,2.791,3.067,3.067,3.728,3.728,4.051,4.051,4.199,4.199,...
            4.450,4.450,4.626,4.626,4.816,4.816,4.960,4.960,5.088,5.088]';
       
% Maturities quoted in table 1
T_mat = (2:10)';
     
%% a) Calibration

% Tenor Structure, T_m = m/2, m = 0,...,19
T = (.0:.5:9.5)';

% For any fixed beta, the system of equations leads to a unique solution
v =  @(beta) Imp_Vol.* sqrt(2*beta ./ (1-exp(-2*beta*T(2:end))).*T(2:end));

% Resulting volatility function for a given set of indices
sigma = @(t,beta) v(beta) .* exp(- beta * (T(2:end) - t));

%% b) MC simulation

% Setup 
M = 20; delta = .5; d_t = 1/12; kappa = .035; K = 1e4;

% Specification I for the diffusion 
lbda_I = @(beta) @(t) sigma(t,beta);

% Specification II for the diffusion
lbda_II = @(beta) @(t) diag(sigma(t,beta));

% T_M - Zero-Coupon Bond Price P(0,T_M) (needed later)
P_0_M = 1/prod(1 + delta * L_0);

% Create a table to store the obtained prices 
Table_b = table(); Table_b.Maturity = T_mat;

Table_b.Quoted = [25.0,77.0,148.5,230.5,325.5,431.5,545.5,664.0,786.0]';

% 3D Matrix containing the cap prices
P = zeros(length(T_mat),2,2);

P(:,:,1) = Cap_Prices(L_0,P_0_M,kappa,K,M,d_t,delta,lbda_I(.07),T,T_mat);

P(:,:,2) = Cap_Prices(L_0,P_0_M,kappa,K,M,d_t,delta,lbda_II(.07),T,T_mat);

s = ["Q_star","Q_M"];

for i = 1:2
    for j = 1:2
        Table_b{:,"Spec_"+num2str(i) + "_" + s(j)} = P(:,j,i);
    end
end

% Display the results
Table_b

%% c) 4 x 6 Swaption

% Reduce the number of MC simulations (too costly otherwise)
K = 1e3;

% Maturity of the swaption and underlying swap respectively
T_mu = 4; T_nu = T_mu + 6;

% Array of beta's
Beta = [-eye(6);(1:6)==(6:-1:1)']*[.4000;.3273;.2545;.1818;.1091;.0364];

% Table containing the approximated prices for Q^* and Q^T_M
Table_c = table(Beta);

% 3D matrix storing the results
Q = zeros(length(Beta),2,2);

init = true;

%===================== ATM 4 x 6-swaption prices =========================%

for b = 1:length(Beta)
        
    if init
        
        % Store the swap rate at the first iteration
        [Q(b,:,1),swp_rate] = Swpt_Price(L_0,P_0_M,T_mu,T_nu,K,M,d_t,delta,...
                                         lbda_I(Beta(b)),T);

        fprintf("Swap Rate: %2.4f\n",swp_rate); init = false;
   
    else
        Q(b,:,1) = Swpt_Price(L_0,P_0_M,T_mu,T_nu,K,M,d_t,delta,...
                              lbda_I(Beta(b)),T,swp_rate);
    end   
        Q(b,:,2) = Swpt_Price(L_0,P_0_M,T_mu,T_nu,K,M,d_t,delta,...
                              lbda_II(Beta(b)),T,swp_rate);                         
    
end
 
% Store the computed prices into a summarising table
for i = 1:2
    for j = 1:2
        Table_c{:,"Spec_"+num2str(i) + "_" + s(j)} = Q(:,j,i);
    end
end

% Display the results
Table_c 

%===================== Plot the swaptions prices =========================%
figure

plot(Beta, Table_c{:,[2,4]},'-' ,"Linewidth",1.5); hold on
plot(Beta, Table_c{:,[3,5]},'--' ,"Linewidth",1.5);

leg = ["Specification I","Specification II"];
legend([leg + ", Q^*",leg + ", Q^{T_M}"],'NumColumns',2)

xlabel("\beta");          ylabel("Swaption price [bps]")
title("Swaption price as a function of \beta")

%% d) General correlation specification

% Use again 10^4 simulations
K = 1e4;

% Focus the computation on the "calibrated" beta, equal to 0.07
Beta = 0.07;

% Array of gamma's
Gamma = [0.1,1,2];

% Correlation matrix as a function of the parameter gamma
rho = @(gamma) exp(- gamma * abs(T(2:end) - T(2:end)'));

% Derive lambda_m(t) = l_m * sigma_m(t) with beta = 0.07
lbda = @(gamma) (@(t) chol(rho(gamma))' .* sigma(t,.07));

% Create a table to store the results
Table_d = table(Beta);

for g = 1:length(Gamma)

    Price_tmp = Swpt_Price(L_0,P_0_M,T_mu,T_nu,K,M,d_t,...
                           delta,lbda(Gamma(g)),T,swp_rate);
                        
    for j = 1:2
        % Store the prices 
        Table_d{:,"Inter_"+num2str(g) + "_" + s(j)} = Price_tmp(j);
    end
                    
end

% Display the results
Table_d


