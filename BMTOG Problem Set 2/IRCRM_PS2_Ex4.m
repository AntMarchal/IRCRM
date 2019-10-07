%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================== Problem Set 2 ============================% 
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all 
clear
warning('off')

%% Exercise 4

%=============================== Part (a) ================================%

T = 5; sigma_ATM = 0.4945; delta = 1; F_4_5 = 0.04;

P_5 = exp(- 0.03 * T);

% ATM Caplet Price 
d = @(i,K,sigma) (log(F_4_5/K) + (-1)^(1+i) * sigma^2 / 2 * (T-1))...
                                            /(sigma * sqrt(T-1));
                                
black_price = @(K,sigma) delta * P_5 * ( F_4_5 * normcdf(d(1,K,sigma))...
                                       - K     * normcdf(d(2,K,sigma)));

Cpl_ATM =  black_price(F_4_5,sigma_ATM)

% Bachelier volatility:
sigma_bach = Cpl_ATM / (P_5 * sqrt(T-1) * normpdf(0))

%=============================== Part (b) ================================%

% Caplet strikes
K = 1e-2 * [1,1.5,2,3,3.5,4,5,6,7,8,9]';

D = (F_4_5 - K) / (sigma_bach * sqrt(T-1));

% Bachelier Prices
Price_bach = delta * P_5 * sigma_bach * sqrt(T-1)...
           * (D .* normcdf(D) + normpdf(D))
       
% Black Implied Volatility
sigma_black = zeros(length(K),1);

for i = 1: length(K)
    f = @(sigma) black_price(K(i),sigma) - Price_bach(i);
    sigma_black(i) = fzero(f,sigma_bach);
    
end
sigma_black
%Implied Black volatility from the original table
IV=1e-2*[64.00,59.65,56.30,51.95,50.50,49.45,47.95,47.00,46.35,45.95,45.70]'
% Plots
plot(K,sigma_black,'o-','linewidth',2)
xlabel('K'); ylabel('Implied black volatility')
hold on
plot(K,IV,'o-','linewidth',2)
hold on
plot([0.01,0.09],[0.4945,0.4945],'-')
legend('bachelier','observe','ATM caplet')

%========================== Part(d) ======================================%

sigma = 0.4945; beta=0.75; delta=1; tau=4; P= exp(- 0.03 * 5);F = 0.04;
K = 1e-2 * [1,1.5,2,3,3.5,4,5,6,7,8,9]';
CplDD=zeros(length(K),1);
for i=1:length(K)
    CplDD(i)=BlackCaplet(delta,P,F/beta,K(i)+F*(1/beta-1),tau,beta*sigma);
end
    

IVBlack=zeros(length(K),1);
for i=1:length(K)
    f=@(sigma) BlackCaplet(delta,P,F,K(i),tau,sigma)-CplDD(i);
    IVBlack(i)=fzero(f,0.49);
end

plot(K,IVBlack,'o-','linewidth',2)
legend('bachelier','observe','ATM caplet','Displaced Diffusion')    



function CplBlack=BlackCaplet(delta,P,F,K,tau,sigma)
d1=(log(F/K)+0.5*sigma^2*tau)/(sigma*sqrt(tau));
d2=(log(F/K)-0.5*sigma^2*tau)/(sigma*sqrt(tau));
CplBlack=delta*P*(F*normcdf(d1)-K*normcdf(d2));
end