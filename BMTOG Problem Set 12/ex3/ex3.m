clear all
%% b) 
global r P0 ;
r = 0.077/0.96; % initial value for r (not specified in problem set)
P0 = p0(1);

delta = 0.5;

Jhandle = @(x) J(x);
JJ = quadgk(Jhandle,0,1);
P1 =  I(1) + (1-delta) * JJ;
fprintf('P0:  %f\n', P0)
fprintf('J:  %f\n', JJ)
fprintf('I:  %f\n',I(1))
fprintf('P1:  %f\n', P1)


%% c)
P0 = p0(5);
Jhandle = @(x) J(x);
JJ = quadgk(Jhandle,0,5);
denom = 0;
for i=1:10
    denom = denom + 0.5 * I(i/2);
end
x = (delta * JJ ) / denom;
fprintf('numer*:  %f\n', delta*JJ)

fprintf('denom*:  %f\n', denom)

fprintf('x*:  %f\n', x)
