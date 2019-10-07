function [res] = hedged_of_s(s)
%HEDGED_OF_S This function reutrns the value
% of the hedged portfolio with a shift s of the yield curve.
global yields hedged_cash_flows
n = length(hedged_cash_flows);
res = hedged_cash_flows.' * exp(- (yields(1:n) + s) .* transpose(1:n));
end

