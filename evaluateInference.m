function [ results ] = evaluateInference( infFunc, synth, opts )
    %evaluateInference A wrapper to evaluate the accuracy of intensity function
    %inference over synthetic data.
    %   
    % Args:
    % * infFunc: a handle to an inference function (the only included one in
    %   this code is runIntensityInference).
    % * synth: a structure containing the synthetic data, known parameters and
    %   known intensity function, constructed by one of the createSyntheticData*
    %   functions. See createSyntheticData for details.
    % * opts: inference options as required by infFunc.
    %
    % Returns:
    % * a structure containing the results of infFunc, plus new fields:
    %   * runtime: the total runtime of the inference
    %   * accuracy: the accuracy of the inference as determined by
    %     evaluateIntensityAccuracy (see for details).
    %
    % Usage:
    % results = evaluateInference(@runIntensityInference, mySyntheticData, myInferenceOpts)
    %
    % See also:
    % runIntensityInference createSyntheticData createSyntheticDataFromEquation
    % evaluateIntensityAccuracy
    
    % (c) Copyright 2013 Thomas A. Lasko
   
    if(isempty(opts.tmax))
        opts.tmax = synth.events(end);
    end
    
    tic;
    results = infFunc(synth.events, opts);
    results.runtime = toc;
    
    results.accuracy = evaluateIntensityAccuracy(results.x, results.intensity, ...
    synth.xIntensity, synth.intensity * synth.params.gammaB / synth.params.gammaA);    
end

