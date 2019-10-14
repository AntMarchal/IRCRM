%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================== Problem Set 4 ============================%
%================================ Exercise 1 =============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format short; warning('off')

%% 0. Setup

% Load the quoted prices
opts = spreadsheetImportOptions("NumVariables", 11);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A22:K352";

% Specify column names and types
opts.VariableNames = ["Date", "CHFConfBonds2y", "CHFConfBonds3y", "CHFConfBonds4y", "CHFConfBonds5y", "CHFConfBonds7y", "CHFConfBonds10y", "CHFConfBonds20y","CHFConfBonds30y", "EURGermanBonds10y", "USDUSBonds10y"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
data = readtable("statmon_E4_M1_M.xls", opts, "UseExcel", false);

% Yield in July 2015
yield = data(331,2:end-2);

% Array of time 
time = [0,2,3,4,5,7,10,20,30]';
T = (0:0.1:30);
yield = [time(2:end),yield{1,:}'];


%% I. Representer Theorem Smoothing Method
% Three different alphas to test
alphas = [0.01 0.1 1];

% Start plot
fig = figure(1);
plot(yield(:,1),yield(:,2),'o')
xlabel('T')
ylabel('Yield Rate'); 
Legend{1} = sprintf('Data');
hold on

% Loop over the three alphas
for s = 1:3
    
    % Scalar Product <h_ih_j>
    hh = zeros(length(time),length(time));
    
    for l = 1:length(time)
        for k = 1:length(time)
            hh(l,k) = time(l)*time(k) + min(time(l),time(k))*time(l)*time(k) ...
                - 1/2*min(time(l),time(k))^2*(time(l)+time(k)) ...
                + 1/3*min(time(l),time(k))^3;
        end
    end

    % Matrix A - Found by creating a linear system of equation from the 
    % Represention Theorem (slide 168 -IRCRM_slides_20190903.pdf )
    A = zeros(9,9);
    % First row: Ti
    for i =1:9
        A(1,i)= 1*time(i);
    end
    % First column: alpha* Ti
    for i = 1:9
        A(i,1)=alphas(s)*time(i);
    end
    % Rest : if i==j, 1+ alpha * <h_i,h_j>, alpha * <h_i,h_j> otherwise
    for i =2:9
        for j =2:9
            if i==j
                A(i,j)= 1+alphas(s)*hh(i,j);
            else
                A(i,j)= alphas(s)*hh(i,j);
            end
        end
    end

    % Array of result
    C = zeros(9,1);
    for i = 2:9
        C(i,1) = alphas(s)*yield(i-1,2)*time(i);
    end
    
    % Array of betas
    Betas = A\C;


    for k = 1:8
        for u = 1:300
            if u>=0
                % smoothed splines
            h(k,u) = time(k+1) + time(k+1)*min(u*0.1,time(k+1)) - 1/2*min(u*0.1,time(k+1))^2;
            else
            h(k,u)= 0;
            end
        end
    end

    % Forward curve function
    F(1)=Betas(1);
    for u =2:300
        F(u) = Betas(1)+Betas(2:9)'*h(:,u-1);
    end
    % yield curve - PAS JUSTE
    for u =1:300
        Y(u) = 1/u*F(1)+Betas(2:9)'*sum(h(:,u),2);
     end
 
    Legend{s+1} = sprintf('alpha = %2.2f',alphas(s))
    plot(T(2:end),F,'-')
end
legend(Legend,'Location','Best'); 

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
[sortedSums, sortOrder] = sort(sum(D, 1), 'Descend')
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
fig(2) = figure
plot(eigenvectors(:,1:3))
legend({'level','slope','curvature'},'Location','Best')
title('PCA Analysis of Monthly Changes')
ylabel('Loading')



