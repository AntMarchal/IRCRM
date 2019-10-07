%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================== Problem Set 2 ============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all 
clear
warning('off')

%% Exercise 1

%=============================== Part (b) ================================%

% Forward curve stored in a vector
Fwd_curve = 1e-2.*[6,8,8,9,9,10,10,11];

% Cash flow dates and time step
T = (1:8)'/4; delta = 1/4;

% Vector of zero-coupon bond prices using the formula derived in the report
fprintf('\nZCB prices (for T = 0,1/4,1/2,...,2):\n')

P = 1 ./ cumprod(1 + delta * Fwd_curve)'

% Computation of the swap rate using the formula from the slides
R_swap = (P(1)-P(end))/(delta * sum(P(2:end)));

fprintf('\nSwap rate: %2.6f\n',R_swap)

%% Exercise 2 

%=============================== Part (a) ================================%

% Nominal, coupon payments and dates 
N = 100; c = 5.5; T_2 = T(2:2:end);

% Vector of cash flows (including the nominal)
cash_flows = c + N * floor(T_2/2);

% Coupon bond price (add the nominal for the maturity date)
p = sum(P(2:2:end) .* cash_flows) ; 

fprintf('\nCoupon bond price: %2.4f\n',p)

%=============================== Part (b) ================================%

% The Yield to Maturity (YTM) of the bond can be obtained by finding 
% the zero of the following function:
g = @(y) sum(exp(-y*T_2) .* cash_flows) - p;

% Choose y_0 = 0 as initial guess:
YTM = fzero(g,0);
fprintf('\nYield to Maturity: %2.4f\n',YTM)

%=============================== Part (c) ================================%

% The expression for the yield curve (YC) in terms of the ZCB prices 
% is derived in the report 

fprintf('\nYield Curve: \n')
YC = - log(P) ./ T 

% We also plot the yield curve
plot(T,YC,'o-','linewidth',2); xticks(T); xlabel('T'); ylabel('y(0,T)')

%=============================== Part (d) ================================%

% We simply use the definitions of each quantity from the lecture

% Macaulay duration
D_Mac = sum(exp(-YTM * T_2) .* T_2 .* cash_flows) / p;

fprintf('\nMacaulay duration: %2.4f\n',D_Mac)

% Duration
D = sum(P(2:2:end) .* T_2 .* cash_flows) / p;

fprintf('\nDuration: %2.4f\n',D)

% Convexity
C = sum(P(2:2:end) .* T_2.^2 .* cash_flows) / p;

fprintf('\nConvexity: %2.4f\n',C)

%=============================== Part (e) ================================%

% Coupon bond price subject to a shift of the yield curve
% (Note that only the yields at T = 1/2,1,3/2,2 contribute to the price)
p_ = @(s)  sum(exp(-(YC(2:2:end)+s) .* T_2) .* cash_flows);

% Plot
figure; ax = axes;

% Choose a set of colors
colors = {[0,0.5,0.8],'k',[0.9,0.5,0],[0.8,0,0]};

fplot(p_,[-0.2,0.2],'linewidth', 2,'color',colors{1}); hold on
xlabel('s'); ylabel('p(s)');

% Add a y tick for the bond price p
yticks([70:10:100,p,110:10:150])
yticklabels({'70','80','90','100','p','110','120','130','140','150'})

%=============================== Part (e) ================================%

style = {'-','--','--'};

for i = 1:3
    
    % Approximation function
    f = @(s) p*(1 - (D + (D_Mac-D) * (i==1)) * s + 1/2 * C * s^2 * (i==3));
    
    % Plot
    fplot(f,[-0.2,0.2],style{i},'linewidth', 2,'color',colors{i+1}); 
    xlabel('s'); ylabel('p(s)')
end

% Highlight the Coupon Bond Price 
plot(0,p,'k.','Markersize',20); y_lim = get(ax,'Ylim'); 
plot([0,0],[y_lim(1),p],'k--','linewidth',0.7)
plot([-0.2,0],[p,p]    ,'k--','linewidth',0.7)

legend('Exact','Order 1 with Macaulay duration','Order 1','Order 2')