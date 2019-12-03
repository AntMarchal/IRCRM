%=========================================================================%
%================= Interest Rate and Credit Risk Models ==================% 
%============================ Problem Set 10 =============================%
%============================== Exercise 1 ===============================%
%=========================================================================%
%======== BRODARD Lionel, MARCHAL Antoine, TISSOT-DAGUETTE Valentin ======%
%======================= OUYANG Tonglin, GIRO Tomas ======================%
%=========================================================================%

clear all; clc; rng(0);

data = readtable('credit.csv', 'delimiter', ';');  % skips the first three rows of data
data = data(:,2:end); % remove first column
T = table2array(data); % convert table to array
X = T(:,[2,3]);
student = T(:,4);
Y = T(:,1);
X_S = X(student==1,:); % _S stands for student
X_R = X(student==0,:); % _R stands for rest ( non student)
Y_S = Y(student==1);
Y_R = Y(student==0);

% In matlab the first class is 1 and the second is 2.
% For us the first class is 1 (default) and the second is 0 (not default).
% For this reason we apply the transformations x:->2-x to our 
% classes to obtain the good order of classes in matlab for the
% mnrfit function.

[B_S,dev_S,stats_S] = mnrfit(X_S, 2-Y_S);
[B_R,dev_R,stats_R] = mnrfit(X_R, 2-Y_R);

disp('latex table of the results of 1a)');
table_res = array2table([B_S, stats_S.se, stats_S.t, stats_S.p,B_R, stats_R.se, stats_R.t, stats_R.p]...
    , 'VariableNames',{'beta_S','sigma_S','z_S', 'p_S', 'beta_R','sigma_R','z_R', 'p_R'}...
    , 'RowNames',{'intercept', 'balance', 'income'});
table_res
table2latex(table_res, 'RegResults.tex');

%[p(B_R,X_R) Y_R]

%% b)
disp('Proba of default of the particular student:');
p(B_S, [500 100])

%% c)
% In matlab the first class is 1 and the second is 2.
% For us the first class is 1 (default) and the second is 0 (not default).
% For this reason we apply the transformations x:->2-x to our 
% classes to obtain the good order of classes in matlab for the
% confusionmat function.
disp('Conf mat of student:');
p_S = p(B_S, X_S)>0.5;
p_R = p(B_R, X_R)>0.5;
cm_S = [sum( (Y_S==1).* (p_S==1)) sum((Y_S==0).*(p_S==1)) ;...
        sum( (Y_S==1).* (p_S==0)) sum((Y_S==0).*(p_S==0)) ];
cm_S = array2table(cm_S, 'VariableNames',{'Real_1', 'Real_0'}, 'RowNames', {'Predicted_1','Predicted_0'})
table2latex(cm_S, 'cm_S.tex');

disp('Conf mat of not student:');
cm_R =  [sum( (Y_R==1).* (p_R==1)) sum((Y_R==0).*(p_R==1)) ;...
         sum( (Y_R==1).* (p_R==0)) sum((Y_R==0).*(p_R==0)) ];
cm_R = array2table(cm_R, 'VariableNames',{'Real_1', 'Real_0'}, 'RowNames', {'Predicted_1','Predicted_0'})
table2latex(cm_R,'cm_R.tex');
