function [L] = MC_Simulation(L_0,M,d_t,delta,lbda,T)

%==========================================================%
% Monte Carlo Simulation - LIBOR Market Model              %
%==========================================================%
% INPUT:                                                   %
% L_0    : Initial LIBOR Forward Curve                     %
% M      : Number of Tenor Dates                           %
% d_t    : Step used in the Euler-Maruyama scheme          %
% delta  : Length of inter tenor dates period              %
% lbda   : Diffusion Function @(t)                         %
% T      : Array of tenor dates                            %
% f      : Payoff function                                 %
%==========================================================%
% OUTPUT :                                                 %
% L      : M x M matrix of simulated LIBOR forward rates   %                                         %
%==========================================================%

%% 0. Setup

% Total number of nodes on the time grid
N = floor(max(T)/d_t);

% Time grid 
T_G = d_t * (0:N);

% Dimension of lambda_m(t)
[~,d] = size(lbda(1));

%% I. Monte Carlo
    
% 3D matrix containing simulated log LIBOR forward rates under Q^* & Q^T_M
H = zeros(M,N+1,2); H(:,1,:) = repmat(log(L_0),1,2);

% White noise
Z = randn(d,N);

for i = 1:N
    
    % Relevant indices to consider (m s.t Tm >= t_i)
    id = find(T > T_G(i));
    
    % Diffusion term (extend the array to get 20 rows)
    lbd = [zeros(1,d);lbda(T_G(i))];

    % Auxiliary array to simplify the notation
    Y = lbd(id,:) .* exp(H(id,i)) ./ (delta^(-1) + exp(H(id,i)));
    
    % Correction term for Q^T_M
    sigma_0_M = sum(Y);

    % Drift, with extra term for Q^T_M (if j = 1)
    alpha = @(j) diag(lbd(id,:) * (cumsum(Y) - sigma_0_M * j)')...
          - 1/2 * sum(lbd(id,:).^2,2);
         
    for j = 1:2   
        % Euler-Maruyama Scheme
        H(id,i+1,j) = H(id,i) + alpha(j-1)*d_t+ lbd(id,:)*Z(:,i)*sqrt(d_t);
    end      
end

% Keep only values at dates T_m, m = 0,...,M-1
L = exp(H(:,ismember(T_G,T),:));

