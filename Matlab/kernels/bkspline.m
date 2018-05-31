function yvals = bspline(xvals, order)
%
%  yvals = bkspline(xvals, order)
%
%  Creates the bspline of order 'order' using the exact 
%  expression from [Unser, Splines: A perfect fit...]. Evaluates 
%  the spline on 'xvals'.
%
%  [Erik Jonsson, 2006]


yvals = zeros(size(xvals));
for k = 0:order+1
  yvals = yvals + nchoosek(order+1, k) * (-1)^k * pospow(xvals-k+(order+1)/2, order);
end
yvals = 1/factorial(order)*yvals;


function y = pospow(x, n)
y = (x>=0) .* (x.^n);

