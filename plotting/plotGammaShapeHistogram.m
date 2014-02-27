function [ ] = plotGammaShapeHistogram( draws, known, color, xLimits )
    %plotGammaShapeHistogram Plots a histogram of the inferred gamma shape
    %parameters.
    %
    % Args:
    % * draws: a vector of parameter values to be displayed.
    % * known: The known correct parameter value.
    % * color: an RGB vector or color name for the histogram fill color.
    % * xLimits: the x-axis limits to be displayed. Optional.
    %
    % Returns nothing.
    %
    % Usage: normally called by plotResults.
    %
    % See also:
    % plotResults
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    setDefault('xLimits', []);
    
    hist(draws, 40);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',color,'EdgeColor',color)
    yLimits = ylim;
    
    if(~isempty(xLimits))
        xlim(xLimits);
    else
        % set xlim to include 0 and known, if provided
        xLimits = xlim;
        xLimits(1) = 0;
        if(~isempty(known) && known > xLimits(2))
            xLimits(2) = known;
        end
        xlim(xLimits);
    end
    
    % Read back into xLimits, just in case matlab did something funny with xlim.
    xLimits = xlim;
    
    if(~isempty(known))
        tickLocations = unique([0, 1, known, xLimits(end)]);
    else
        tickLocations = get(gca, 'XTick');
    end
    
    set(gca, 'FontSize', 8, 'XTick', tickLocations, 'LineWidth', 0.0001);
    text(0.5, -0.25, 'Shape Param', 'Units', 'normalized', ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'middle')
    box off;
    set(gca, 'Color', 'none');
    set(gca, 'YTick', []);
    set(gca, 'YColor', [.7 .7 .7]);
    
    line([1 1], ylim, 'Color', [.3 .3 .3]);
    line(xLimits, [0 0], 'Color', 'black');
    
    if(~isempty(known))
        line([known known], ylim, 'Color', 'red');
    end

end

