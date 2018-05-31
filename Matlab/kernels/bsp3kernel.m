function y = bsp3kernel(x)
%
%  y = bsp3kernel(x)
%
%  The 3rd order b-spline kernel, evaluated on the values in 'x'. 
%  [Erik Jonsson, 2006]

x = abs(x);

y = (x<1) .* (2/3 - x.^2 + x.^3/2) + ...
    (x>=1 & x < 2) .* (2-x).^3/6;

    