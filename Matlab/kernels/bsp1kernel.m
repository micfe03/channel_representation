function y = bsp1kernel(x)
%
%  y = bsp1kernel(x)
%
%  The first order B-spline kernel, evaluated at the points in x.
%  [Erik Jonsson, 2006]

y = (abs(x) < 1) .* (1-abs(x));
