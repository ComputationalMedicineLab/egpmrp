function [yAxisMax] = plotResults( results, synth, yAxisMax, fig, fillColor, uncertainty )
    %PLOTRESULTS plots the inference results in a nice way.
    %
    % Computes the median and 95%CI (if requested by uncertainty=true) from
    % results, then plots using plotIntensity with an inset by
    % plotGammaShapeHistogram.
    %
    % Args:
    % * results: the structure produced by runIntensityInference (see for
    %   details) or evaluateInference.
    % * synth: the structure produced by one of the createSyntheticData*
    %   functions (see for details).
    % * yAxisMax: the maximum value for the y axis. Omit or set to [] for auto.
    % * fig: a figure handle. Omit or set to [] to create a new figure.
    % * fillColor: an rgb vector representing the color to fill space between
    %   the inferred intensity function and the x-axis. Omit or set to [] for no fill.
    % * uncertainty: a boolean indicating whether to draw 95% CIs. 
    %
    % Returns:
    % * The yAxisMax that was used in the plot, for setting future plots to
    %   match.
    %
    % See also:
    % plotIntensity plotGammaShapeHistogram
    
    % (c) Copyright 2013 Thomas A. Lasko

    setDefault('yAxisMax', []);
    setDefault('fillColor', []);
    setDefault('uncertainty', false);
    
    if(~isempty(synth))
        known = [synth.xIntensity, synth.intensity / synth.params.gammaA];
        knownGammaA = synth.params.gammaA;
    else
        known = [];
        knownGammaA = [];
    end
    
    if(uncertainty)
       p = [50, 2.5, 97.5];
       iqi = prctile(results.intensity, p, 1);
       
       mid = iqi(1,:);
       low = iqi(2,:);
       hi = iqi(3,:);
   else
       p = 50;
       iqi = prctile(results.intensity, p, 1);
       
       mid = iqi(1,:);
       low = [];
       hi = [];
    end     

   if(~isempty(fillColor))
       histColor = darken(fillColor);
   else
       histColor = [.6 .6 1];
   end
   
    plotIntensity(results.x, results.events, mid, hi, low, fig, known, yAxisMax, fillColor);
    mainAxes = gca;
    
    xlabel('Time');
    ylabel('Intensity');
    ylimits = ylim;
    yAxisMax = ylimits(2);

    if(~isempty(results.param))
        axes('Position',[.78 .78 .2 .2])
        %color = [.6 .6 1];
        plotGammaShapeHistogram([results.param.gammaA], knownGammaA, histColor);
    end
    
    % Set axes back to main axes so subsequent post-plot-manipulation happens to
    % the main plot.
    set(gcf, 'CurrentAxes', mainAxes);
    
    hold off;
end

function [lightColor] = lighten(color) %#ok<DEFNU>
    lightColor = (color + [1, 1, 1]) / 2;
end

function [darkColor] = darken(color)
    darkColor = (color + [0, 0, 0]) / 2;
end
