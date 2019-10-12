function [p] = compute_zcoupon_price(t)
%COMPUTE_ZCOUPON_PRICE 
global b0 b1 b2 a
% We compute the zero coupon bond price
p = t.*b0+b1.*(1-exp(-a.*t))+b2./a*(1-exp(-a.*t)+a.*t.*exp(-a.*t));
end

