function [ events ] = sampleEventsFromIntensity( xIntensity, intensity, gammaA, gammaB )
    %sampleEventsFromIntensity Summary of this function goes here
    %   Detailed explanation goes here
    delta = xIntensity(2) - xIntensity(1);
    warp = intensityToWarp(intensity, delta);
    events = sampleWarpedGamma(xIntensity, warp, gammaA, gammaB);
end

function [ events ] = sampleWarpedGamma(xWarp, warp, gammaA, gammaB) 
    %sampleWarpedGamma Sample from a warped gamma distribution
    %  Samples first from a homogeneous gamma distribution, then
    %  warps the event times into the non-homogeneous data space.
    
    rate = gammaA * gammaB;
    
    % Create inter-event intervals for twice the expected number of events
    numSamples = ceil(2 * (warp(end) - warp(1)) / rate);
    intervals = gamrnd(gammaA, gammaB, numSamples, 1);
    
    % First interval is from the incomplete homogeneous gamma pdf
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
