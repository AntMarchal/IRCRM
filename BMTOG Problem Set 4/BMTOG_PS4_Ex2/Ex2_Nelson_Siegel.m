%=========================================================================%
%============== Interest Rate and Credit Risk Models =====================% 
%============================== Problem Set 4 ============================%
%================================ Exercise 2 =============================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

close all; clear; clc; format short; warning('off')

%% 0. Setup

% Load the quoted prices
Quotes = readtable("Bootstrap_data.xls",'Range','C5:C22');

Quotes.Properties.VariableNames = {'Dates'};

% Load the exact discount ratesobtained in the previous assignment
% using the pseudo-inverse method
Table = readtable("Discount_curves.csv");

id = ismember(Table.Dates,Quotes.Dates);

DC = Table.DC_pseudo_inverse(id);

% Spot Date
t_0 = datetime('03-Oct-2012');

% Vector ot time to maturities of the quoted prices
T_quoted = delta(t_0,Quotes.Dates);

% Yield curve
fprintf('\nYield curve for the quoted dates only:\n');
YC = - log(DC)./ delta(t_0,Quotes.Dates)

%% I. Nelson-Siegel Smoothing Method

% Parameters
a = 2e-2 * (3:5); 

% Neslon-Siegel (NS) yield curve 
NS_YC = @(a,beta,T) beta(1) + (beta(2) + beta(3))...
         .* (1 - exp(- a * T))./(a * T) - beta(3) .* exp(- a * T);
     
% Forward curve function
f = @(a,beta,T) beta(1) + (beta(2) + beta(3)*(a * T)).* exp(- a * T);
     
% Function returning the squared error
Squared_error = @(a,beta) norm(NS_YC(a,beta,T_quoted) - YC,2)^2;
     
% First guess for beta_i, i = 0,1,2
beta_0 = 1/2 * ones(3,1);

Squared_errors = zeros(3,1); beta_opt = ones(length(a),3);

% Plot
figure; Legend = cell(1,length(a));

% Styles and colors
style = {'-','--','-.'}; color = {[0,0.8,1],[0.9,0,0.5],'k'};

grid = linspace(0,max(T_quoted),1000);

for i = 1:length(a)

    fprintf('\n=================== a = %2.2f ===================\n',a(i))
    % Optimal parameters obtained by minimizing the squared error
    beta_opt(i,:) = fminsearch(@(beta) Squared_error(a(i),beta),beta_0);

    fprintf('\nOptimal Parameters:\n'); beta_opt(i,:)

    Squared_errors(i) = Squared_error(a(i),beta_opt(i,:));
    
    fprintf('\nSquared error: %e\n',Squared_errors(i));     
    
    plot(grid,f(a(i),beta_opt(i,:),grid),style{i},...
        'Color',color{i},'Linewidth',1.5); hold on
    
    Legend{i} = sprintf('a = %2.2f',a(i));
        
end

xlabel('Time to Maturity'); xlim([0,max(T_quoted)]); 
ylabel('Forward Rate [%]'); 

legend(Legend,'Location','Best'); title('Forward Curves')

%% II. Compute and plot the discount curve for the tenor dates

% Time to Maturity for all the tenor dates
Tenor_Dates = delta(t_0, Table.Dates);

% Neslon-Siegel (NS) yield curve 
NS_YC = @(a,beta,T) beta(1)  + (beta(2) + beta(3))...
         .* (1 - exp(- a * T))./(a * T) - beta(3) .* exp(- a * T);

% Formula to pass from yield curve to discount curve
YC_to_DC = @(a,beta,T) exp(- NS_YC(a,beta,T) .* T);

% Plot 
figure

for i = 1:length(a)
    
    % Column name
    name_tmp = "DC_smoothing_a_" + num2str(100*a(i)) + "_percent";
    
    Table{:,name_tmp} = YC_to_DC(a(i),beta_opt(i,:),Tenor_Dates);
    
   
    plot(Tenor_Dates,Table{:,name_tmp},style{i},'Color',color{i},...
         'Markersize',3,'Linewidth',1.5); hold on


    
    Legend{i} = sprintf('a = %2.2f',a(i));
end

xlabel('Time to Maturity'); xlim([0,max(Tenor_Dates)]); 
ylabel('Discount Rate [%]'); 

legend(Legend,'Location','Best'); title('Discount Curves')

Table
