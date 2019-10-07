classdef Bond
    %BOND This object represents a bond
    properties
        cash_flows
        price
        maturities
        duration
        convexity
    end
    
    methods
        function obj = Bond(cash_flows)
            %BOND Constructs an instance of this class. 
            %   This function is called when the new bond is initialized
            %   and computes the price, duration and convexity.
            global zcbonds
            obj.cash_flows = cash_flows;
            obj.maturities = transpose(1:length(obj.cash_flows));
            % The price is the sum of the cash flows times the z-coupon
            % bonds
            obj.price = obj.cash_flows.'*zcbonds(1:length(obj.cash_flows));
            % The duration is the inverse of the price times the sum of the
            % cash flows multiplied by the z-coupon bonds and the
            % maturities.
            obj.duration = 1/obj.price * (obj.cash_flows.'*(obj.maturities.*zcbonds(1:length(obj.cash_flows)) ));
            % The convexity is the inverse of the price times the sum of the
            % cash flows multiplied by the z-coupon bonds and the
            % maturities squared.
            obj.convexity = 1/obj.price * (obj.cash_flows.'*((obj.maturities.^2).*zcbonds(1:length(obj.cash_flows)) ));
        end
            
        
       
        
       
    end
end

