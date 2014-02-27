function [ warp ] = intensityToWarp( intensity, samplingPeriod )
    %IntensityToWarp Integrates an intensity function to a warp function.
    %   Uses trapezoid integration, assuming that intensity is measured at
    %   regular times spaced by samplingPeriod.
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    warp = zeros(size(intensity));
    
    % Trapezoid integration, assuming regular sampling of x
    warp(2:end) = cumsum((intensity(1:end-1) + intensity(2:end)) * samplingPeriod /2);
    
end

