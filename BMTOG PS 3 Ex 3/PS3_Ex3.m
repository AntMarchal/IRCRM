%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================== Problem Set 3 ============================%
%================================ Exercise 3 =============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format short; warning('off')

%=========================================================================%
%============================== a) Bootstrap  ============================%
%=========================================================================%

%% 0. Setup

% Load the data 
Data = readtable("Bootstrap_data.xls",'Range','C5:E22');

% Change the column names
Data.Properties.VariableNames = {'Dates' 'Quotes' 'Source'};

% Transform the percentages into decimal values
Data.Quotes = Data.Quotes/100

% Table that will contain the discount curve
Table = readtable("Bootstrap_data.xls",'Range','H5:H44');

% Change the column name
Table.Properties.VariableNames = {'Dates'};

% Add a column for the discount curve (DC) obtained with bootstrap
Table.DC_bootstrap = zeros(size(Table));

% Spot Date
t_0 = datetime('03-Oct-2012');

% Indices at which a quoted rate is available for a given instrument
id = @(W) find(ismember(Table.Dates,W));

%% I. LIBOR

% Indices corresponding to LIBOR 
l = strcmp(Data.Source,'LIBOR');

% Maturity Dates for LIBOR
S = Data.Dates(l);

% Quoted LIBOR rates
L_S = Data.Quotes(l);

% Use the relationship between ZCB prices and interest rates
Table.DC_bootstrap(id(S)) = (1 + delta(t_0,S) .* L_S ).^(-1);

%% II. Futures

% Indices corresponding to Futures
f = strcmp(Data.Source,'Futures');

% Dates of the quoted futures prices Maturities
T = Data.Dates(f); 

% Add the first reset date (T_1)
T = [T(1) - calmonths(3) + caldays(1); T];

% Array containing P(t_0,T_i) 
P_T = zeros(length(T),1);

% Linear interpolation to obtain L(t_0,T_1) (see linear_interpol.m)
% (note that S_3 < T_1 < S_4)
L_T_1 = linear_interpol(S(3:4),L_S(3:4),T(1));

% Deduce P(t_0,T_1)
P_T(1) = (1 + delta(t_0,T(1)) * L_T_1)^(-1);

% Forward rates (using futures as proxy)
F = 1 - Data.Quotes(f);

% Use the command diff() to store all the durations delta(T_i,T_i+1)
Delta_T = days(diff(T))/360;

% Use the iterative formula to compute P(t_0,T_i) for all i > 2
P_T(2:end) = P_T(1)./ cumprod(1 + Delta_T.* F);

% Add P_T to the Table and obtain the whole discount curve
Table.DC_bootstrap(id(T)) = P_T;

%% III. Swaps

% Indices corresponding to Swaps
s = strcmp(Data.Source,'Swap');

% Dates of the quoted swap rate 
U_quoted = Data.Dates(s);

% Array of coupon dates for the swaps
U =  [Table.Dates(id(U_quoted(1))) - calyears(1); 
      Table.Dates(id(U_quoted(1)):end)         ];
           
% Array for the discount curve P(t_0,U_k) and vector of swap rates
P_U = zeros(length(U),1); R_swap = P_U;

% Shifted indices at which a quoted swap rate is available
id_U = id(U_quoted) - (length(S) + length(T));

% Add the quoted swap rates in R_swap
R_swap(id_U) = Data.Quotes(s);

%======================= Derivation of P(t_0,U_2) ========================%

% Linear interpolation to obtain L(t_0,U_1) (note that T_3 < U_1 < T_4)

% First derive the LIBOR rates L(t_0,T_i), i = 3,4
L_T_2_3 = (1./Table.DC_bootstrap(id(T(3:4))) -1)./delta(t_0,T(3:4));

L_U_1 = linear_interpol(T(3:4) , L_T_2_3 , U(1));

% Deduce P(t_0,U_1) 
P_U(1) = (1 + delta(t_0,U(1)) * L_U_1)^(-1);

% Compute R_swap(t_0,U_1) (required for the pseudo-inverse approach)
R_swap(1) = (1 / P_U(1) -1)/delta(t_0,U(1));

% Derive P(t_0,U_2) using the inverted swap rate formula:
P_U(2) = (1 - R_swap(2) * delta(t_0 ,U(1)) * P_U(1))...
       / (1 + R_swap(2) * delta(U(1),U(2)));

%============ Linear Interpolation of the missing swap rates =============%

for j = 1:length(id_U)-1
    
    % indices between index id(k) and id(k+1) and boundaries
    id_tmp = (id_U(j)+1):(id_U(j+1)-1); id_b = [id_U(j),id_U(j+1)];

    % Linear interpolation of the missing swap rates 
    R_swap(id_tmp) = linear_interpol(U(id_b),R_swap(id_b),U(id_tmp));

end

%=================== Compute P(t_0,U_k) iteratively ======================%

% Store the increments delta(U_n-1,U_n), n > 0 (U_0 = t_0)
Delta_U = [delta(t_0,U(1)); days(diff(U))/360];

for n = 3:length(U)
    
    % Use the iterative formula
    P_U(n) = (1 - R_swap(n) * dot(Delta_U(1:(n-1)),P_U(1:(n-1))))...
           / (1 + R_swap(n) * Delta_U(n));
end

% Add P_U to the Table and obtain the whole discount curve
Table.DC_bootstrap(id(U)) = P_U;

%=========================================================================%
%=================== b) Pseudo-inverse on increments =====================%
%=========================================================================%

%% 0. Construction of Matrix of Cash Flows C and array of prices p

% Array of prices (simply add ones when the quote is a LIBOR or swap rate)
p = l + s; 

n = length(p); N = length(Table.Dates); C = zeros(n,N); 

%=========================== LIBOR cash flows ============================%
 
 for k = 1:length(S)
     C(k,ismember(Table.Dates,S(k))) = 1 + delta(t_0,S(k)).* L_S(k);
 end
 
%=========================== Futures cash flows ==========================%
 
 for k = 1:length(T)-1     
     C(k + length(S),id(T(k:(k+1)))) = [-1, 1 + delta(T(k),T(k+1)) * F(k)];
 end
     
%======================= Swap cash flows (U_0 = t_0) =====================%
 
% First row index in C to add the swap cash flows
shift = length(S) + length(T) - 1;

 for k = 1:length(U_quoted)
     
     % Column indices corresponding to the cash flows
     id_tmp = id(U(U <= U_quoted(k)));

     % Add the cash flows (note that c_k = 1 + delta * K) 
     C(k + shift,id_tmp) = Delta_U(1:id_U(k)).* R_swap(id_U(k))...
                         + (id_tmp == id_tmp(end));
 end

%% I. Pseudo-inverse solution 

% Array of time increments delta(x_i,x_{i-1}) = delta(t_i,t_{i+1})
delta_ = days(diff(Table.Dates))/360;

% Set up the matrices M, W, e_1 and A  
M = eye(N) - diag(ones(N-1,1),-1);

W = diag(1./sqrt([delta(t_0,S(1)); delta_]));

e_1 = (1:length(M) == 1)';

A = C * (M \ inv(W));

% Estimation of the vector of weighted increments
Delta_star = A' * ((A * A') \( p - C * (M \ e_1)));

% Resulting discount curve
Table.DC_pseudo_inverse = M \ (W \ Delta_star + e_1)

%=========================================================================%
%================================ Plots  =================================%
%=========================================================================%

% Array of Time to Maturity
TTM = delta(t_0,[t_0;Table.Dates]);

%============================ Discount Curves ============================%
 
figure

plot(TTM,[1;Table.DC_bootstrap],'o-','Linewidth',2,'MarkerSize',6); hold on
plot(TTM,[1;Table.DC_pseudo_inverse],'k:','Linewidth',2); 

legend('Bootstrap','Pseudo-inverse')
xlabel('Time to Maturity'); ylabel('P(t_o,T)');
axis([0,TTM(end) + 1,0.4,1.1])

title('Discount Curves')

%====================== Forward and LIBOR curves =========================%
 
% One should add F(t_0,t_0,S_1) = L(t_0,S_1)

% Function computing the forward curve for a given discount curve
Fwd_curve = @(DC) 100.* [L_S(1); (1/DC(1)-1)/delta(t_0,Table.Dates(1));
                                 (DC(1:end-1)./DC(2:end)- 1)./ delta_];
                    
% Function computing the LIBOR curve for a given discount curve
LIB_curve = @(DC) 100 .*[L_S(1); (1./DC - 1)./ TTM(2:end)];

dc = {'DC_bootstrap','DC_pseudo_inverse'}; style = {'o-','o:'};

figure

for i = 1:2
    
plot([TTM,TTM],[Fwd_curve(Table{:,dc{i}}),LIB_curve(Table{:,dc{i}})],...
               style{i},'MarkerSize',4,'Linewidth',2); hold on

fprintf('\nForward rate between 10/2041 and 10/2042 for %s: %2.4f %%\n',...
         dc{i},Fwd_curve(Table{:,dc{i}})'* (TTM == TTM(end)))

end

xlim([0,TTM(end) + 1]); xlabel('Time to Maturity'); ylabel('Rate [%]')

legend('Forward rate - Bootstrap','LIBOR rate    - Bootstrap',...
       'Forward rate - Pseudo-inverse','LIBOR rate    - Pseudo-inverse',...
       'Location','SouthEast')

title('Implied LIBOR and Forward Curves')