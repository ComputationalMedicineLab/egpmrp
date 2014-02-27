function [ ] = xRug( x, tickSize, color )
    %xRug Draws rug plot with ticks at (x, 0).
    % Intended to be added to an existing plot.
    % Args:
    % * x: a vector of locations for ticks
    % * tickSize: the size of ticks in the units of ylim. Defaults to 
    %   .01 * (ymax - ymin).
    % * color: an rgb vector or color name for ticks. Defaults to 'red'.
    %
    % Returns nothing.
    %
    % Usage:
    % plot(x, y);
    % xRug(x);
    %
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    setDefault('tickSize', .01);
    setDefault('color', 'red');
    
    yAxisLimits = ylim;
    tickLen = tickSize * diff(yAxisLimits);
    if(size(x, 2) == 1)
        x = x';
    end
    xd = repmat(x, 3, 1);
    xd(1,:) = nan;
    yd = repmat([yAxisLimits(1), yAxisLimits(1), yAxisLimits(1)+tickLen], 1, length(x));
    plot(xd(:), yd(:), '-', 'Color', color);
end

