classdef Bond
    %BOND This object represents a bond
    properties
        cash_flows
        price
        b0sens
        b1sens
        b2sens
    end
    
    methods
        function obj = Bond(cash_flows);
            %BOND Constructs an instance of this class. 
            %   This function is called when the new bond is initialized
            %   and computes the price and sensitivity to each beta
            global zcbonds maturities a
            obj.cash_flows = cash_flows;
            obj.price = cash_flows*zcbonds.';
            obj.b0sens = (cash_flows.*(-maturities))*zcbonds.';
            obj.b1sens = (cash_flows.*(-1-exp(-a.*maturities)))*zcbonds.';
            obj.b2sens = (cash_flows.*(-(1/a.*(1-exp(-a.*maturities)+a.*maturities.*exp(-a.*maturities)))))*zcbonds.';
        end

    end
end

