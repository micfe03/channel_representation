function yvals = bspkernel_diff(xvals)
%
%  yvals = bspkernel_diff(xvals)
%
%  

% First option
% yvals = bkspline(xvals+0.5, 1) - bkspline(xvals-0.5, 1);

% Second option
absx = abs(xvals);
yvals = (absx < 1/2) .* (-2*xvals) + ...
        (xvals >= 1/2 & xvals <= 3/2) .* (xvals-1.5) + ...
        (xvals <= -1/2 & xvals >= -3/2) .* (xvals+1.5);

% It's easy to verify that they give the same results
