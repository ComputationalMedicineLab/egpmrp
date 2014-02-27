function [ warp, intensity ] = GpToWarp(gpX, gpY)
    %GpToWarp Converts a GP (a log intensity curve) to a warping curve.
    
    % (c) Copyright 2013 Thomas A. Lasko
   
    intensity = GpToIntensity(gpY);
    
    % infer the spacing between grid points, assumed to be constant
    delta = gpX(2) - gpX(1);
    
    warp = intensityToWarp(intensity, delta);
    
end

