function [ res ] = p( B, X )
%P computes the probability
M = size(X,1);
X = [repmat(1, M, 1) X];
x = X*B;
res = exp(x) ./ (1+exp(x));
end

