function [ opts ] = intensityInferenceDefaults()
    %intensityInferenceDefaults Returns a struct containing default parameters for intensity inference.
    %
    % Args:
    % None
    %
    % Returns:
    % * opts: a struct with default initialization. See runIntensityInference for
    % details.
    %
    % Usage:
    % opts = intensityInferenceDefaults()
    % runIntensityInference(data, opts);
    %
    % See also:
    % runIntensityInference

    % (c) Copyright 2013 Thomas A. Lasko
    
    opts.numIntensityPoints = 200;
    opts.numBurnin = 1000;
    opts.numIter = 1000;
    opts.thetaUpdateProb = 0.5;
    opts.silent = false;
    opts.fig = [];
    opts.tmin = 0;
    opts.tmax = []; % default use last event as tmax
    opts.thresholdNumEvents = 3;

end

