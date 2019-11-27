function [Prices] = Cap_Prices(L_0,P_0_M,kappa,K,M,d_t,delta,lbda,T,T_mat)

%==========================================================%
% Interest Cap Pricing  - LIBOR Market Model               %
%==========================================================%
% INPUT:                                                   %
% L_0    : Initial LIBOR Forward Curve                     %
% P_0_M  : Initial T_M ZCB Price (P(0,T_M))                %
% kappa  : Cap Rate                                        %
% K      : Number of Monte Carlo Simulations               %
% M      : Number of Tenor Dates                           %
% d_t    : Step used in the Euler-Maruyama scheme          %
% delta  : Length of inter tenor dates period              %
% lbda   : Diffusion Function @(t)                         %
% T      : Array of tenor dates                            %
% T_mat  : Array of cap maturities                         %
%==========================================================%

% Payoff function
f = @(L) delta * max(diag(L) - kappa,0);

% Array of simulated prices for different maturities
Pi = zeros(M-1,K,2);

for k = 1:K
    
    % Monte Carlo Simulation of forward LIBOR rates
    L = MC_Simulation(L_0,M,d_t,delta,lbda,T);

    %======================== Under Q^* ========================% 
    
    % Implied Money Market Account B^*(T_m)
    B_star = cumprod(1 + delta * diag((L(:,:,1))));

    % Discount the cash flows w.r.t. the money market account
    Pi(:,k,1) = cumsum(f(L(2:end,2:end,1))./B_star(2:end));

    %====================== Under Q^{T_M} ======================%
    
    % Implied Zero Coupon Prices P(T_m,T_M), m = 2,...,M
    P_m_M = [1./prod(1 + delta*tril(L(3:end,3:end,2)),1),1]';

    % Use Lemma 37.2 (Change of Numéraire)   
    Pi(:,k,2) = cumsum(f(L(2:end,2:end,2))./P_m_M) * P_0_M ;  
end

% Resulting prices (in basis points)
Prices = 1e4 * permute(mean(Pi,2),[1 3 2]);

% Keep only cap prices for maturities in 2,...,10 years
Prices  =  Prices(ismember(T(2:end) + delta,T_mat),:);