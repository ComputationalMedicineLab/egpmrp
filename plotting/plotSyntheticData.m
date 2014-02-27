function [] = plotSyntheticData(data, fig)
    %PlotSyntheticData A convenience function for plotting synthetic data.
    % plotSyntheticData(data, fig) plots the synthetic intensity function with
    % overlaid events contained in data. Uses plotIntensity.
    %
    % Args:
    % * data: a structure of synthetic data, such as generated from 
    %   createSyntheticData (see for details).
    % * fig: a figure handle. Defaults to a new figure.
    %
    % Returns nothing.
    %
    % See also:
    % createSyntheticData createSyntheticDataFromEquation plotIntensity
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    setDefault('fig', @()(figure()));
    plotIntensity([], data.events, [], [], [], fig, [data.xIntensity, data.intensity]);
end

