function [ intensity ] = GpToIntensity(gp)
    %GpToIntensity Converts a GP to an intensity curve.
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    intensity = exp(gp) + 1e-12; % add small correction to avoid zero intensity

end

