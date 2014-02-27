function [ results ] = inferIntensityAndParameters( model, opts)
    %inferIntensityAndParameters draws intensity function and model
    % parameters from the posterior density given the events and priors.
    %   This functin runs a specified number of MCMC iterations inferring both the
    % intensity function and the model parameters.
    
    % Args:
    % * model: a structure containing model information necessary for the MCMC
    %   inference, constructed by createModel (see for details).
    %   
    % * opts: a structure of inference options, with the following fields
    %   * tmin, tmax: the beginning and ending of inference time
    %     interval. Must contain all event times, but may extend beyond them.
    %     Setting opts.tmax = [] will set tmax to be the last event time.
    %   * numIntensityPoints: The number of evenly-spaced grid points on
    %     which to define the intensity function. Run time is cubic in this
    %     parameter, but it must be large enough to capture the desired
    %     variations in the intensity function.
    %   * thresholdNumEvents: Must have more than this number of events to
    %     perform the inference. Otherwise returns a constant function of
    %     appropriate intensity.
    %   * numBurnin, numIter: Number of burn-in and inference
    %     iterations.
    %   * silent: Turns off periodic progress updates and plotting.
    %   * fig: If set to a figure handle, plots the intensity function and
    %     prints parameter values at each iteration and prints para. Useful for
    %     understanding what's going on, but slows down the inference quite a
    %     bit. Overridden by opts.silent=true.
    %   * thetaUpdateProb: The probability that a given iteration will update
    %     the GP parameters as well as the intensity function, instead of the
    %     intensity function alone. thetaUpdateProb near 1 causes almost every
    %     iteration to update parameters (faster convergence per iteration), and
    %     thetaUpdateProb near 0 causes almost no iteration to update parameters
    %     (faster iteration). The default of 0.5 is a reasonable starting choice.
    %
    % Returns:
    % * a struct containing information from each of the m = opts.numIter
    %   iterations. Each component of the struct is a matrix or vector with the
    %   information from iteration i contained in row i:
    %   * intensity: an matrix of width n = opts.numIntensityPoints 
    %     containing draws of the intensity function values corresponding to the 
    %     locations at model.xIntensity
    %   * ll: a vector containing the posterior log likelihood of the
    %     intensity function and model parameters given the events.
    %   * param: a vector containing structs of the model parameters: 
    %     * covHyp: the log length scale and log magnitude variance
    %       parameters for the GP covariance function.
    %     * gammaA: the gamma shape parameter
    %
    % Usage: Generally called by runIntensityInference
    % results = inferIntensityAndParameters(model, opts);
    %
    % See also:
    % intensityInferenceDefaults createModel
    
    % (c) Copyright 2013 Thomas A. Lasko
    
    
    x = model.events;
    n = size(x, 1);
    
    if (n <= opts.thresholdNumEvents)
        constIntensity = n / (model.xIntensity(end) - model.xIntensity(1));
        results.intensity = constIntensity * ones(size(model.xIntensity))';
        results.param = [];
        results.ll = [];
        return;
    end
    
    numIntensityPoints = size(model.xIntensity, 1);
    
    % initialization
    intensityDraws = zeros(numIntensityPoints, opts.numIter);
    llDraws = -Inf * ones(opts.numIter, 1);
    
    paramDraws(1:opts.numIter) = struct('covHyp', [0; 0], ...
        'gammaA', 0);
    param = getInitialParams(model);
    
    KFun = makeKFun(model.xIntensity, model.covFunc);
    cholSigma = computeCholSigma(KFun, param.covHyp);
    
    llFun = makeLogLikFun(model.events, model.xIntensity, param.gammaA);
    
    gp = cholSigma' * randn(numIntensityPoints, 1);
    ll = llFun(gp);
    
    % run inference
    for i=(1-opts.numBurnin):(opts.numIter)
        try
            if(rand(1) < opts.thetaUpdateProb)
                [gp, param.covHyp, cholSigma] = updateGpAndParams(gp, param.covHyp, ...
                    model.gpParamPriorFun, KFun, llFun);
                model.hypCov = param;
            else
                [gp, ll] = updateGp(gp, ll, llFun, cholSigma);
            end
            
            [ll, llFun, param.gammaA] = updateGammaA(gp, model.events, model.xIntensity, ...
                param.gammaA, model.gammaAPriorFun, llFun);
            
            if(~opts.silent && (mod(i, 100) == 0 || ~isempty(opts.fig)))
                fmt = '%d: ell: %5.5g, sf: %5.5g gamma a: %2.2g, ll %g.\n';
                fprintf(fmt, i, exp(param.covHyp(1)), exp(param.covHyp(2)), ...
                    param.gammaA, ll);
            end
            
            if(~opts.silent && ~isempty(opts.fig))
                figure(opts.fig);
                clf;
                plot(model.xIntensity, gp, 'ob', 'MarkerSize', 1.5);
                hold on;
                plot(x, pchip(model.xIntensity, gp, x), '.r', 'MarkerSize', 10);
                ylim([-10, 5]);
                xlim([floor(model.xIntensity(1)), ceil(model.xIntensity(end))]);
                set(gca,'LooseInset',get(gca,'TightInset'));
                drawnow;
            end
            
            if(i > 0)
                intensityDraws(:,i) = GpToIntensity(gp) / param.gammaA;
                llDraws(i) = ll/n; % Avg log likelihood per point in x.
                paramDraws(i) = param;
            end
        catch exception
            fprintf('Error in iteration %d: %s %s\n', i, exception.identifier, exception.message);
        end
    end
    results.intensity = intensityDraws';
    results.ll = llDraws;
    results.param = paramDraws;
end


function [ llf ] = makeLogLikFun( events, xIntensity, gammaA )
    %MAKELOGLIKFUN Returns a function that computes the posterior log likelihood
    %of an intensity function and the gammaA parameter given the events.
    
    x = events;
    n = size(x, 1) - 1;
    K = n * (-gammaln(gammaA));
    
    function [ll] = f(gp)
        % Returns the posterior log likelihood of the log intensity function gp
        
        [w, intensity] = GpToWarp(xIntensity, gp);
        
        if(isinf(w(end)))
            ll = -Inf;
            return
        end
        
        % Project the data x into the warped space
        
        % pchip maintains the monotonicity of the warp function. spline does not.
        wx = pchip(xIntensity, w, x);
        tail = pchip(xIntensity, w, [x(end); xIntensity(end)]);
        head = pchip(xIntensity, w, [xIntensity(1); x(1)]);
        
        % Compute intensities at the events (plus realmin to avoid zero intensity).
        ix = pchip(xIntensity, intensity, x) + realmin;
        
        dx = diff(wx) + realmin; % avoid zero differences
        dt = diff(tail);
        dh = diff(head);
        
        
        assert(dt >= 0 && dh >= 0, 'Bug detected - dt: %f, dh: %f', dt, dh);
        if(dt == 0)
            llTail = 0;
        else
            llTail = -dt; % Poisson approximation for tractability
        end
        
        if(dh == 0)
            llHead = 0;
        else
            llHead = -dh; % Poisson approximation for tractability
        end
        ll = K - sum(dx) + (gammaA - 1) * sum(log(dx))  + sum(log(ix)) + llTail + llHead;
    end
    llf = @f;
end

function [ ll, llFunUpdate, gammaAUpdate ] = updateGammaA( gp, events, xIntensity, gammaA, ...
        gammaAPriorFun, llFun)
    
    logstep = 0.1;
    
    % Default is no update
    llFunUpdate = llFun;
    gammaAUpdate = gammaA;
    
    ll = llFun(gp);
    logA1 = log(gammaA);
    plogA1 = gammaAPriorFun(logA1);
    
    % lognormal update
    logAStep = logstep*randn(1);
    logA2 = logA1 + logAStep;
    plogA2 = gammaAPriorFun(logA2);
    
    gammaA2 = exp(logA2);
    llFun2 = makeLogLikFun(events, xIntensity, gammaA2);
    
    ll2 = llFun2(gp);
    alpha = ll2 - ll + logAStep + plogA2 - plogA1;
    
    if(log(rand) < alpha)
        llFunUpdate = llFun2;
        ll = ll2;
        gammaAUpdate = gammaA2;
    end
end


function [fUpdate, covHypUpdate, cholSigmaUpdate] = ...
        updateGpAndParams( f, covHyp, covHypPriorFun, KFun, llFun )
    %UPDATEGPANDPARAMS Updates f and covHyp using surrogate data slice sampling.
    
    aux = 1; % Stdev of noise to add to f.
    sliceWidth = 10;
    
    [covHypUpdate, fUpdate, ~, cholSigmaUpdate] = update_theta_aux_surr(...
        covHyp, f, llFun, KFun, aux, covHypPriorFun, sliceWidth);
end

function [ fUpdate, llUpdate ] = updateGp( f, ll, llFun, cholSigma )
    %UPDATEGP Updates f using elliptical slice sampling
    
    [fUpdate, llUpdate] = elliptical_slice(f, cholSigma, llFun, ll, 0);
end

function [ KFun ] = makeKFun(xLocs, covfunc)
    %MAKEKFUN Create function to compute covariance matrix at xLocs given covariance parameters.
    cf = ensureCell(covfunc);
    x = xLocs;
    eps = 1e-6 * eye(size(x, 1)); % Small diagonal offset for numeric stability
    
    function [K] = f(params)
        % Assumes params is a column vector of the appropriate log(params)
        % for the covariance function in the model.
        K = feval(cf{:}, params, x);
        if(K(1,1) > 1000)
            fprintf('Large K: %f\n.', K(1,1));
        end
        K = K + eps;
    end
    KFun = @f;
end

function [params] = getInitialParams(model)
    initLogEll = unifrnd(log(model.ellMin), log(model.ellMax));
    initLogSf = log(unifrnd(model.sfMin, model.sfMax));
    
    params.covHyp = [initLogEll; initLogSf];
    params.gammaA = unifrnd(model.gammaAMin, model.gammaAMax);
end

function [ cs ] = computeCholSigma(KFun, covHyp)
    K = KFun(covHyp);
    kc = 1e-6 * eye(size(K,1)); % tiny correction for numeric stability.
    cs = chol(K + kc);
end