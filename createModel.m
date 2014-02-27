function [ model ] = createModel(events, tmin, tmax, numIntensityPoints)
    %createModel creates a struct containing model parameters.
    %   Modifying this file is the way to specify alternate priors or GP 
    % covariance functions.
    %
    % Returned model structure contains the following components:
    %   * events: a vector of event times to be modeled. Just a copy of the
    %     input events vector, for convenience in passing to other functions.
    %   * xIntensity: a vector of numIntensityPoints locations/times at which to 
    %     calculate the intensity function, equally spaced from tmin to tmax.
    %   * sfMin, sfMax: min and max allowed values for the GP magnitude
    %     variance parameter.
    %   * ellMin, ellMax: min and max allowed values for the GP length scale.
    %   * gpParamPriorFun: a one-argument handle to a function giving the prior
    %     probability of a vector of GP parameter values. Currently truncated
    %     exponential for the length scale, uniform for magnitude.
    %   * covFunc: the GP covariance function.
    %   * gammaAMin, gammaAMax, min and max allowed values for the gamma shape
    %     parameter.
    %   * gammaAPriorFun: a one-argument handle to a function giving the prior
    %     probability of a given gamma shape parameter. Currently uniform over
    %     [gammaAMin, gammaAMax].
    
    % (c) Copyright 2013 Thomas A. Lasko

    setDefault('numIntensityPoints', 200);
    
    model.events = events;

    % Intensity function sampling locations
    model.xIntensity = linspace(tmin, tmax, numIntensityPoints)';
    
    % Gaussian process parameter priors
    model.sfMin = 0.01;
    model.sfMax = 7;
    model.ellMin = 2 * tmax / numIntensityPoints;
    model.ellMax = 2 * tmax;
    
    model.gpParamPriorFun = makeSEPriorFun(model.ellMin, model.ellMax, ...
        model.sfMin, model.sfMax);

    % Gaussian process covariance function 
    model.covFunc = @covSEiso;
        
    % Gamma process parameter priors and functions
    model.gammaAMin = 0.2;
    model.gammaAMax = 25;
    model.gammaAPriorFun = makeUniformPriorFun(model.gammaAMin, model.gammaAMax);    
    
end


% building blocks for parameter hyperprior functions
function [ lpFun ] = makeSEPriorFun(ellMin, ellMax, sfMin, sfMax)
    %%MAKESEPRIORFUN returns a prior function on [log(ell), log(sf)].
    
    % prior used for sinewave synthetic data, see paper for reasoning.
    % ellPriorFun = makeTruncatedPriorFun(makeLogGaussianPriorFun(log(.3), .3), ellMin, ellMax);
    
    % prior for other synthetic data, see paper for reasoning.
    ellPriorFun = makeTruncatedPriorFun(makeExpPriorFun(10), ellMin, ellMax);
    % other priors are possible, such as the unused 1/x prior below.
    
    sfPriorFun = makeUniformPriorFun(sfMin, sfMax);
    
    function [logProb] = prior(logParams)
        logProb = ellPriorFun(logParams(1)) + sfPriorFun(logParams(2));
    end
    
    lpFun = @prior;
end

function [ priorFun ] = makeTruncatedPriorFun( f, xMin, xMax )
%MAKETRUNCATEDFUN Truncates prior function f to return -Inf outside [xMin, xMax].
    unif = makeUniformPriorFun(xMin, xMax);
    
    function [logProb] = prior(logX)
        logProb = unif(logX) + f(logX);
    end
    priorFun = @prior;
end

function [ priorFun ] = makeExpPriorFun( mean )
    %MAKEEXPPRIORFUN Returns function f(x) that computes log prob x under
    %exponential prior.
    a = 1/mean;
    loga = log(a);
    
    function [logProb] = prior(logX)
        x = exp(logX);
        logProb = loga - x * a + log(abs(logX)) + logX;
    end
    
    priorFun = @prior;
end

function [priorFun] = makeUniformPriorFun(xMin, xMax)
    % MAKEUNIFORMPRIORFUN Creates function f(x) computing the prior of x under
    % uniform probability.
    logMin = log(xMin);
    logMax = log(xMax);
    
    function [logProb] = prior(logX)
        if(logMin <= logX && logX <= logMax)
            logProb = 0;
        else
            logProb = -Inf;
        end
    end
    
    priorFun = @prior;
end

function [ priorFun ] = makeLogGaussianPriorFun( logMean, sd ) %#ok<DEFNU>
    %MAKELOGGAUSSIANPRIORFUN Returns function f(x) that computes log prob x under
    %log gaussian prior.
       
    function [logProb] = prior(logX)
        logProb = log(normpdf(logX, logMean, sd));
    end
    
    priorFun = @prior;
end

function [ priorFun ] = makeInversePriorFun() %#ok<DEFNU>
    %MAKEPOWERPRIORFUN Returns function f(x) that computes log prob x under
    %1/x prior (modulo the constant).
       
    function [logProb] = prior(~)
        logProb = 0;
    end
    
    priorFun = @prior;
end

