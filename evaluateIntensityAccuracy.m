function [ accuracy ] = evaluateIntensityAccuracy( xIntensity, inferredIntensity, xKnown, knownIntensity )
    %evaluateIntensityAccuracy evaluates the accuracy of intensity function
    %inference on synthetic data, by comparing with the known generating
    %intensity function.
    %
    % Args:
    % * xIntensity: a vector of times of the intensity function values. 
    % * inferredIntensity: a matrix in which each row corresponds to an inferred
    %   intensity function, measured at times xIntensity. Each row must have the
    %   same width as xIntensity.
    % * xKnown: a vector containing the known intensity function times.
    % * knownIntensity: a vector containing the values of the known intensity
    %   function, measured at times xKnown. Must be the same length as xKnown.
    %
    % Returns:
    % * a structure containing fields:
    % * mse: The mean squared error of all inferred intensity functions compared
    %   to knownIntensity, calculated at the points xIntensity);
    % * logP: The log probability of the known intensity function given the
    %   inferred intensity function density, calculated by summing at each point
    %   xIntensity. The intensity function density at each point is calculated by
    %   kernel smoothing (in other words, it is not assumed to be Gaussian).
    %
    % Usage:
    %   Usually called by evaluateInference, but directly callable if desired:
    %   accuracy = evaluateIntensityAccuracy(xIntensity, inferredIntensity, xKnown, knownIntensity);
    %
    % See also:
    % evaluateInference
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    interpKnownIntensity = spline(xKnown, knownIntensity, xIntensity);
    mse = 0;
    logP = 0;
    n = size(inferredIntensity, 2);
    
    for i = 1:n
        % normal kernel with theoretically optimal bw
        pd = fitdist(inferredIntensity(:, i), 'Kernel');
        
        logP = logP + max(-800, log10(pdf(pd, interpKnownIntensity(i))));
        mse = mse + (interpKnownIntensity(i) - icdf(pd, 0.5))^2;
    end
    accuracy.mse = mse / n;
    accuracy.logP = logP;
end

