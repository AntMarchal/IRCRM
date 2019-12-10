%% Interest Rate and Credit Risk Models
% Problem Set 11
% OUYANG Tonglin, BRODARD Lionel, MARCHAL Antoine
% TISSOT-DAGUETTE Valentin, GIRO Tomas
clc
clear all
%% a)
global L r
B = 50;
V = 100;
r = 0.05;
T = 5;
sigma = 0.25;
dp = ( log(V/B) + ( r + sigma^2/2 ) * T ) / (sigma * sqrt(T) );
dm = ( log(V/B) + ( r - sigma^2/2 ) * T ) / (sigma * sqrt(T) );
S_0 = V * normcdf(dp) - exp(-r* T ) * normcdf(dm) * B;
disp('The value of the equity is: ');
disp(S_0);

B_0 = exp(-r*T) * B + V * normcdf(-dp) - exp(-r * T) * normcdf(-dm) * B;
disp('The value of the debt is: ');
disp(B_0);

%% b)
for L=[0.3 0.6 0.9]
     N = 100;
    [T,sigma] = meshgrid(linspace(0,10,N),linspace(0,0.5,N) );
    C = c(T, sigma);
    figure
    surf(T,sigma, C);
    title(strcat('L=',num2str(L)));
    xlabel('T');
    ylabel('volatility');
    zlabel('c');
end

