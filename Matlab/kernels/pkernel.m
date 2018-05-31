function y = pkernel(x)
%
%  y = pkernel(x)
%
%  Evaluate the p-channel kernels on all points in x. Returns a matrix y
%  with twice as many rows as x, since the p-channel representation uses
%  two kernels. The histograms and offset components are interleaved.
%
%  [Erik Jonsson, 2006]


% Histogram component
yh = double(x>=-0.5 & x<0.5);

% Offset component
yo = yh .* x;

% Interleave the two
y = zeros(size(x,1)*2, size(x,2));
y(1:2:end,:) = yh;
y(2:2:end,:) = yo;
