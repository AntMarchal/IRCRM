function [Price,kappa]=Swpt_Price(L_0,P_0_M,T_mu,T_nu,K,M,...
                                  d_t,delta,lbda,T,kappa)

%==========================================================%
% Swaption Pricing - LIBOR Market Model                    %
%==========================================================%
% INPUT:                                                   %
% L_0    : Initial LIBOR Forward Curve                     %
% P_0_M  : Initial T_M ZCB Price (P(0,T_M))                %
% T_mu   : Maturity of the swaption                        %
% T_nu   : Maturity of the underlying swap                 %           
% K      : Number of Monte Carlo Simulations               %
% M      : Number of Tenor Dates                           %
% d_t    : Step used in the Euler-Maruyama scheme          %
% delta  : Length of inter tenor dates period              %
% lbda   : Diffusion function @(t)                         %
% T      : Array of tenor dates                            %
% kappa  : Strike rate                                     %  
%==========================================================%

% Coupon dates of the underlying swap (/!\ annual coupon)
T_c = (T_mu + 1):T_nu;

% Indices corresponding to the coupon dates 
id = find(ismember(T + delta,T_c));

% Indices corresponding to T_mu and T_nu respectively
mu = find(T + delta == T_mu); nu = find(T + delta == T_nu);

% If no strike is given, assume that the swaption is ATM
if nargin < 11
    
    % Initial discount curve P(0,T_m), m = 1,...,20
    P = 1./cumprod(1 + delta * L_0);

    % Swap rate (s.t. the swaption is ATM)
    kappa = (P(mu) - P(nu))/ sum(P(id));
end

% Array of simulated swaption prices
Pi = zeros(1,K);

for k = 1:K
    
% Monte Carlo Simulation of forward LIBOR rates
L = MC_Simulation(L_0,M,d_t,delta,lbda,T);

    for j = 1:2    
        
    % Implied Zero Coupon Prices P(T_mu,T_m), m = mu + 1,...,nu = M
    P = 1./cumprod(1 + delta * L(mu+1:nu,mu,j));
    
    % Payoff (/!\ annual coupon)
    f = max(dot(P(id - mu), L(id,mu,j) - kappa),0);
    
    %======================== Under Q^* ========================% 
    if j == 1
            
        % Implied Money Market Account B^*(T_m)
        B_star_mu = prod(1 + delta * diag((L(1:mu,1:mu,j))));
    
        % Discount the cash flow using the money market account
        Pi(:,k,j) = f / B_star_mu;
        
    %====================== Under Q^{T_M} ======================%
    else
        
        % Use Lemma 37.2 (Change of Numéraire)      
        Pi(:,k,j) = P_0_M * f / P(end);
    end

    end

end

% Resulting prices (in basis points)
Price = 1e4 * permute(mean(Pi,2),[1 3 2]); 
