function [ data ] = createSyntheticDataFromEquation( params )
    %createSyntheticDataFromEquation Draws synthetic data from a
    %function-modulated gamma process.
    %   createSyntheticDataFromEquation( params ) draws events from a modulated
    % gamma process using the warping method, with params.equation as the
    % intensity function. The number of sampled points is determined by the
    % intensity function. The events and the details of the model that generated
    % them are returned. The model details can be used for checking the accuracy
    % of downstream inference methods.
    %
    % Args:
    % * params: a struct of model parameters, containing
    %   * tmin, tmax: the min and max times for the samples.
    %   * equation: a one-argument function handle returning a vector of
    %      intensities for a given vector of locations.
    %   * gammaA: the gamma process shape parameter, must be nonnegative.
    %   * gammaB: the gamma process scale parameter, must be positive.
    %
    % Returns a struct containing the synthetic data, with fields:
    % * events: a vector of synthetic event times
    % * intensity: a vector of intensity values from the intensity function that
    %   generated the events.
    % * xIntensity: a vector of times of the intensity values.
    % * warp: a vector of warp values (cumulative intensity values),
    %    computed using smooth integration of intensity vs. xIntensity.
    % * params: a copy of the input params that generated the events.
    %
    % Usage:
    %  params = syntheticDefaults();
    %  params.equation = @(x) 5 * sin(x.^2) + 6;
    %  data = createSyntheticDataFromEquation(params);
    %
    % See also:
    % syntheticDefaults createSyntheticData
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    numPts = 1000;
    
    data.xIntensity = linspace(params.tmin, params.tmax, numPts)';
    data.intensity  = params.equation(data.xIntensity);
    %     delta = data.xIntensity(2) - data.xIntensity(1);
    %     data.warp = intensityToWarp(data.intensity, delta);
    %     data.events = sampleWarpedGamma(data.xIntensity, data.warp, params.gammaA, params.gammaB);
    data.events = sampleEventsFromIntensity(data.xIntensity, data.intensity, ...
        params.gammaA, params.gammaB);
    data.params = params;
end

