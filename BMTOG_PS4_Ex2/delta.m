function d = delta(t_1,t_2)

%==================================================================%
% Function returning the elapsed time in years between two dates   %
% Convention: actual/360                                           %
%==================================================================%

d = caldays(between(t_1,t_2,'Days'))/360;