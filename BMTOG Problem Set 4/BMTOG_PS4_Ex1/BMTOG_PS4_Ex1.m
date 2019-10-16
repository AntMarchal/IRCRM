%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================== Problem Set 4 ============================%
%================================ Exercise 1 =============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format short; warning('off')

%==================================== (a) ================================%

%% 0. Setup

% Load the quoted prices
opts = spreadsheetImportOptions("NumVariables", 11);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A22:K352";

% Specify column names and types
opts.VariableNames = ["Date", "CHFConfBonds2y", "CHFConfBonds3y",...
                      "CHFConfBonds4y", "CHFConfBonds5y", "CHFConfBonds7y",...
                      "CHFConfBonds10y", "CHFConfBonds20y","CHFConfBonds30y",....
                      "EURGermanBonds10y", "USDUSBonds10y"];
opts.VariableTypes = ["datetime", "double", "double", "double",...
    "double", "double", "double", "double", "double", "double", "double"];

% Import the data
data = readtable("statmon_E4_M1_M.xls", opts, "UseExcel", false);

% Yield in July 2015 
Y = data{end,2:end-2}';

% Array of time to maturities
T = [2,3,4,5,7,10,20,30]'; N_T = length(T);

%% I. Representer Theorem Smoothing Method

% Three different alphas to test
alpha = [0.01 0.1 1]; N_alpha = length(alpha);

% Array of optimal betas
Beta = zeros(length(Y) + 1, N_alpha); 

% Vectorize the linear system of equations on slide 168:
% Find beta in R^{n+1} s.t. A * beta = C 
% (please refer to the report for the definition of A and C)

% Matrix containing the scalar products <h_i,h_j>_H
H = (T * T') .* (1 + min(T,T')) - 1/2 * min(T,T').^2 .* (T + T')...
  + 1/3 * min(T,T').^3;

C = [0;Y.*T];

% Function handle for the quadratic splines
h = @(u,T) T + u .* (T - u/2) .* (u <= T) + 1/2 * T.^2 .* ((u > T) );

%% Derive the betas and yield curve for each alpha

figure

for i = 1:N_alpha

% Construct the matrix A 
A = [[0;T],[T';(H + eye(N_T)/alpha(i))]];

% Solve for beta
Beta(:,i) = A\C;

% Array of integrals of h_i(u) between 0 and some time s > 0
h_int = @(s,T) (T * s + 1/2 * T * s^2 - 1/6 * s^3 ) .* (s <= T)...
      + (s * (T + T.^2 / 2)- 1/6 * T.^3) .* (s >  T);
  
% Resulting yield curve [%]
y = @(s) dot(Beta(:,i), [1;1/s * h_int(s,T)]);

% Plot
fplot(y, [0,T(end)],'Linewidth',1.2); hold on
title('Yield Curves')

Legend{i} = "\alpha =" + num2str(alpha(i));

end

% Add the quoted yields 
plot(T,Y,'d','LineWidth',1.7); 
xlabel('Time to Maturity')  ; ylabel('Yield Rate [%]'); 

Legend{N_alpha + 1 } = sprintf('Data');
legend(Legend,'Location','Best'); 

%% Instantaneous forward curve

% Anonymous function for the forward curve
f = @(u,beta) dot(beta, [1;h(u,T)]);

figure

for i = 1: N_alpha
    
fplot(@(u) f(u,Beta(:,i)), [0,T(end)],'Linewidth',1.2); hold on

title('Instantaneous Forward Curves')

Legend{i} = "\alpha =" + num2str(alpha(i));

end

xlabel('Time to Maturity')  ; ylabel('Instantaneous Forward Rate [%]'); 
legend(Legend,'Location','Best'); ylim([-1,1.4])

%==================================== (b) ================================%

%% II. Monthly Changes
% Number of Months
t_0 = datetime('August-2005');
t_1 = datetime('July-2015');
nMonths=calmonths(between(t_0,t_1,'Month'));

% Monthly Yields used for this exercises
monthlyYield = table2array(data(212:end,2:9));

% Compute monthly change
for i =1:8
    for t = 1:nMonths
        x(t,i)=monthlyYield(t+1,i)-monthlyYield(t,i);
    end
end

% Perform PCA
% Empirical Covariance Matrix
covarianceMatrix = cov(x);
sprintf('Empirical Covariance Matrix:')
display(covarianceMatrix)

% Eigenvalues and Eigenvectors
[V,D]=eig(covarianceMatrix);
[sortedSums, sortOrder] = sort(sum(D, 1), 'Descend');
eigenvectors=V(:,sortOrder); % Eigenvectors
eigenvalues=D(:,sortOrder);% Eigenvalues
display(eigenvectors)
display(eigenvalues)

% Explained Variance
dataInPrincipalComponentSpace = x*eigenvectors;
varianceDataPCS = var(dataInPrincipalComponentSpace)';
explained = 100*varianceDataPCS./sum(varianceDataPCS);
display(explained)

% Plot of the PCA
figure
plot(eigenvectors(:,1:3),'Linewidth',1.2)
legend({'level','slope','curvature'},'Location','Best')
title('PCA Analysis of Monthly Changes')
ylabel('Loading')
