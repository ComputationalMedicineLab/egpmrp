function [ results ] = runIntensityInference( events, opts )
    %runIntensityInference Top-level convenience function to infer model parameters and
    %intensity function from events.
    %
    %   This is the main entry point for inference over real data. To get
    % alternate model behavior such as alternate priors or covariance functions 
    % (as opposed to simply alternate inference options) you should modify
    % createModel or its helper functions.
    %
    % Args:
    % * events: a vector of event times. Times must be represented as scalar 
    %   numeric values, rather than any kind of time/date structure. (For
    %   example, use 10.25 to mean 6:00 AM on day 10.)
    % * opts: a structure of inference options, see inferIntensityAndParameters
    %    for details. Defaults are produced by intensityInferenceDefaults().
    % Returns:
    % * a struct containing all of the fields from inferIntensityAndParameters, 
    %   plus:
    %    * x: a vector of intensity function location/time grid points.
    %    * events: the vector of input event times, for convenience.
    %
    % Usage:
    %  events = getMyEventsSomehow(...); 
    %  opts = intensityInferenceDefaults();
    %  opts.thetaUpdateProb = 0.8; % modify opts as desired
    %  results = runIntensityInference(events, opts);
    %
    % See also:
    % intensityInferenceDefaults inferIntensityAndParameters
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    model = createModel(events, opts.tmin, opts.tmax, opts.numIntensityPoints);
    
    results = inferIntensityAndParameters(model, opts);
    results.x = model.xIntensity;
    results.events = events;
end