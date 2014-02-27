function [] = plotIntensity(xIntensity, events, yMid, yHi, yLow, fig, known, yAxisMax, fillColor)
    %plotIntensity Plots known and inferred intensities and events in a standard way.
    %  This function is used by more convenient, higher level plotting routines.
    %
    % Args:
    % * xIntensity: a vector of times of intensity values in the other args.
    % * events: a vector of original event times.
    % * yMid: a vector of the median inferred intensity values.
    % * yHi: a vector of the upper confidence interval (such as the 97.5th
    %   percentile) of inferred intensity values. Optional.
    % * yLow: a vector of the lower confidence interval (such as the 2.5th
    %   percentile) of inferred intensity values. Optional.
    % * fig: a handle to an existing figure in which to plot. Optional.
    % * known: a vector of the known intensity values. Optional.
    % * yAxisMax: The y axis maximum. Optional.
    % * fillColor: RGB vector of the color to fill between the inferred
    % intensity and the x axis. Optional.
    %
    % Returns nothing.

    % See also:
    % plotResults plotSyntheticData
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    setDefault('fig', @() figure());
    setDefault('xIntensity', []);
    setDefault('yMid', []);
    setDefault('yHi', []);
    setDefault('yLow', []);
    setDefault('known', []);
    setDefault('yAxisMax', []);
    setDefault('fillColor', []);
    
    uncertaintyColor = [0.8 0.8 1];
    lineColor = [1 1 1];
    edgeColor = uncertaintyColor;
    
    if(~isempty(fillColor))
        uncertaintyColor = lighten(fillColor);
        edgeColor = fillColor;
        lineColor = fillColor;
    end
    
    numPlotPoints = 2000;
    if(~isempty(xIntensity))
        xPlot = linspace(xIntensity(1), xIntensity(end), numPlotPoints);
    elseif(~isempty(known))
        xPlot = linspace(known(1,1), known(end,1), numPlotPoints);    
    end
    
    figure(fig);
    hold on;
    if(~isempty(yHi) && ~isempty(yLow))
        yHi = interp1(xIntensity, yHi, xPlot, 'spline');
        yLow = interp1(xIntensity, yLow, xPlot, 'spline');
        if(~isempty(fillColor))
            fill([xPlot, xPlot(end), xPlot(1)], [yHi, 0, 0],  fillColor);
        end
        fill([xPlot, flipdim(xPlot,2)], [yHi, flipdim(yLow, 2)], ...
            uncertaintyColor, 'EdgeColor', edgeColor);
    end
    
    set(gca,'LooseInset',get(gca,'TightInset'));
    xlabel('Time (days)');
    ylabel('Intensity (events/day)');
    
    
    if(~isempty(yMid))
        yMid = interp1(xIntensity, yMid, xPlot, 'spline');
        if(~isempty(fillColor))
            fill([xPlot, xPlot(end), xPlot(1)], [yMid, 0, 0],  fillColor, 'EdgeColor', darken(fillColor));
        else
            plot(xPlot, yMid, 'color', lineColor, 'LineWidth', 1);
        end
    end
    
    if(~isempty(known))
        plot(known(:,1), known(:,2), '-r', 'LineWidth', 1);
    end
        
    if(~isempty(events))
        yLimits = ylim;
        if(isempty(yAxisMax))
            yAxisMax = yLimits(2);
        end
        yAxisMin = -0.03 * yAxisMax;
        ylim([yAxisMin, yAxisMax]);
        h = line(xlim, [0 0], 'Color', 'black'); %#ok<NASGU>
        xRug(events, .02, 'red');
        
        % Move the zero baseline behind the intensity curve
%        uistack(h, 'bottom');
    end
    set(gca, 'TickLength', [0.005, 0.005]);
    hold off;
end

function [lightColor] = lighten(color)
    lightColor = (color + [1, 1, 1]) / 2;
end

function [darkColor] = darken(color)
    darkColor = (color + [0, 0, 0]) / 2;
end

