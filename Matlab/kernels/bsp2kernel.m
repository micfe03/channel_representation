function y = bsp2kernel(x)
%
%  y = bsp2kernel(x)
%
%  The second order B-spline kernel, evaluated at the points in x.
%  [Erik Jonsson, 2006]

y = (abs(x) < 1/2) .* (3/4 - abs(x).^2) + ...
    (abs(x) >= 1/2 & abs(x) <= 3/2) .* (3/2 - abs(x)).^2/2;
