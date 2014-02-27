function [ events ] = sampleWarpedGamma(xWarp, warp, gammaA, gammaB) 
    %sampleWarpedGamma Sample from a warped gamma distribution
    %  The function first samples from a homogeneous gamma distribution, then
    %  warps the event times into the real data space.
    % Args:
    % * xWarp: a vector of times corresponding to warp function values.
    % * warp: a vector of warp function values.
    % * gammaA: the gamma shape parameter.
    % * gammaB: the gamma scale parameter (always = 1 in the paper).
    %
    % Returns:
    % * a vector of event times.
    %
    % Usage: Called by createSyntheticData* functions.
    %
    % See also:
    % createSyntheticData createSyntheticDataFromEquation
    
    % (c) Copyright 2013 Thomas A. Lasko
   
    rate = gammaA * gammaB;
    
    % Create twice the expected number of intervals
    numSamples = ceil(2 * (warp(end) - warp(1)) / rate);
    intervals = gamrnd(gammaA, gammaB, numSamples, 1);
    
    % First interval is from the incomplete gamma pdf
    intervals(1) = gammaincinv(rand(1), gammaA, 'upper') * gammaB;
    
    % Compute events in the warped space from intervals
    wEvents = cumsum(intervals);
    
    % Add more events if necessary
    numExtra = ceil(0.2 * numSamples);
    while(wEvents(end) < warp(end))
        wEvents = [wEvents; wEvents(end) + cumsum(gamrnd(gammaA, gammaB, numExtra, 1))]; %#ok<AGROW>
    end
    
    % Trim the excess
    wEvents(wEvents > warp(end)) = [];
    
    % Unwarp to data space
    events = pchip(warp, xWarp, wEvents);
end


