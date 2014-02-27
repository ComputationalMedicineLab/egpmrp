function [ params ] = syntheticDefaults( )
    %syntheticDefaults Returns a struct containing default parameters for 
    % modulated gamma process event synthesis. 
    %
    % Args:
    %  None
    % 
    % Returns:
    % * params: a struct with default initialization. See createSyntheticData
    %   and createSyntheticDataFromEquation for field details.
    %
    % Usage:
    %  params = syntheticDefaults();
    %  data = createSyntheticData(params);
    %
    % See also:
    % createSyntheticData createSyntheticDataFRomEquation
    
    % (c) Copyright 2013 Thomas A. Lasko

    params.sf = 3.5;
    params.ell = 40;
    params.mean = -6;
    params.skew = 0;
    params.gammaA = 3;
    params.gammaB = 1;
    params.tmax = 1000;
    params.tmin = 0;
    params.equation = @(x) 1;
end

