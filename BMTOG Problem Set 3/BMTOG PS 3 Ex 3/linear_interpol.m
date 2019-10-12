function interpolated_value = linear_interpol(x,y,w)

%===================================================%
% Linear Interpolation of w s.t. x_1 < w < x_2      % 
%                                                   %
% x = [x_1,x_2]: interpolation nodes                %
% y = [y_1,y_2]: images of x_1,x_2                  %             
%                                                   %    
%===================================================%

q = delta(w,x(2))/delta(x(1),x(2));

interpolated_value = [q 1-q] *  y;
