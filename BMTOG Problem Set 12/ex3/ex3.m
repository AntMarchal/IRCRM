clear all
%% b) 
global r ; % r is the initial value for r
r = 0.077/0.96; % initial value for r (not specified in problem set)

delta = 0.5;

Jhandle = @(x) J(x); % we define a handle for the function J
JJ = quadgk(Jhandle,0,1); % we compute the C* from the report
P1 =  I(1) + (1-delta) * JJ; % we compute the price p1
fprintf('P1:  %f\n', P1)


%% c)
Jhandle = @(x) J(x);
JJ = quadgk(Jhandle,0,5);
denom = 0;  % we compute the denominator of the formula on page 677
for i=1:10
    denom = denom + 0.5 * I(i/2); 
end
x = (delta * JJ ) / denom;


fprintf('x*:  %f\n', x)
