function [ data ] = createSyntheticDataFromParams( params )
    %createSyntheticData Draws synthetic data from a Gaussian-process-modulated gamma process.
    %   Creates data by first drawing a squared-exponential Gaussian process
    % (GP) from the prior specified by params, then drawing data points from a
    % modulated gamma process using the warping method, with the GP as the log
    % intensity function. The number of sampled points is determined by the
    % intensity function. The events and the details of the model that generated
    % them are returned. The model details can be used for checking the accuracy
    % of downstream inference methods.
    %
    % Args:
    % * params: a structure with data creation parameters:
    %   * tmin, tmax: the min and max times for the samples.
    %   * ell: the GP length scale
    %   * sf: the GP magnitude (variance) scale
    %   * mean: the constant GP mean 
    %   * gammaA: the gamma process shape parameter, gammaA >= 0
    %   * gammaB: the gamma process scale parameter, gammaB > 0
    %   * skew: if nonzero, adds skew to the GP. A negative skew that is
    %     relatively small compared to sf produces a more realistic
    %     approximation of medical event data. Skewing is an extra step here
    %     that is not replicated in the inference. But it seems to pose no
    %     problem for inference, at least at the values for skew that I've
    %     tried.
    %
    % Returns a struct with fields:
    % * events: a vector of synthetic event times
    % * intensity: a vector of values for the intensity function that generated 
    %   the events. 
    % * xIntensity: a vector of times of the intensity values.
    % * params: a copy of the input params that generated the events.
    %   
    % Usage:
    %  params = syntheticDefaults();
    %  data = createSyntheticData(params);
    %
    % See also:
    % syntheticDefaults createSyntheticDataFromEquation
    
    % (c) Copyright 2013 Thomas A. Lasko

numPts = max(1000, 4 * (params.tmax - params.tmin) / params.ell);

data.xIntensity = linspace(params.tmin, params.tmax, numPts)';
data.intensity = createIntensityFunction(data.xIntensity, params.ell, params.sf, ...
    params.mean, params.skew);
data.events = sampleEventsFromIntensity(data.xIntensity, data.intensity, params.gammaA, params.gammaB);
data.params = params;

end

function [intensity] = createIntensityFunction(xIntensity, ell, sf, mean, skew)
    gp = drawGp(xIntensity, ell, sf, mean);

    % Skew the gp if desired. 
    f = 2 .* gp .* normcdf(skew * gp / sf, 0, 1);
    
    intensity = GpToIntensity(f);
end
    

function [f] = drawGp(x, ell, sf, mn)
    mu = meanConst(mn, x);
    K = covSEiso([log(ell); log(sf)], x);
    kc = 1e-6 * eye(size(x,1));  % tiny correction for numeric stability.
    f = chol(K + kc)' * randn(size(x)) + mu;
end

