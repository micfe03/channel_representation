function yvals = rectkernel(xvals)
%
%  yvals = rectkernel(xvals)
%
%  Evaluates a rectangular box function on the values in 'xvals'.
%  This can be used to create regular histograms in the same code 
%  framework as general channel encodings.
%
%  [Erik Jonsson, 2006]

yvals = (xvals > -0.5) & (xvals <= 0.5);
clss = class(xvals);
yvals = eval([clss '(yvals)']);

