% Top-level script to test the creation of synthetic events using a warped gamma
% process and to test inference from those events. Used to generate results for
% the paper:

% Lasko TA. Efficient Inference of Gaussian Process Modulated Renewal Processes 
% with Application to Medical Event Data. UAI 2014 (and ArXiv:1402.4732 [stat.ML]).

% (c) Copyright 2013 Thomas A. Lasko

%% Create events from gp
%rng(2)
params = syntheticDefaults();

params.sf = 1;
params.ell = 70;
params.mean = -3.0;
params.skew = -0.5;
params.gammaA = 0.5;
params.gammaB = 1;
params.tmax = 3000;
params.tmin = 0;

synth = createSyntheticDataFromParams(params);
fig = figure(1);
clf;
plotSyntheticData(synth, fig);
fprintf('%d events\n', numel(synth.events));

%% Or, create events from equation
params.gammaA = 3;
params.gammaB = 1; % should remain fixed at 1
params.tmin = 0;
params.c = params.gammaA;

params.equation = @(x) params.c * (5 * sin(x.^2) + 6);
params.tmax = 5;
eqLabel = 'sin';

% OR
% params.equation = @(x) params.c * interp1([0, 25, 50, 75, 100], [2, 3, 1, 2.5, 3], x);
% params.tmax = 100;
% eqLabel = 'piecewise-linear';
 
% OR
% params.equation = @(x)  params.c * (2 * exp(-x/5) + exp(-((x - 25)/10).^2)); 
% eqLabel = 'exp';
%params.tmax = 50;

synth = createSyntheticDataFromEquation(params);
fig = figure(1);
clf;
plotSyntheticData(synth, fig);
fprintf('%d events\n', numel(synth.events));


%% Test recovery of known parameters and intensity
trainFig = []; %figure(1);
resultsFig = figure(3);
inferenceOpts = intensityInferenceDefaults();
inferenceOpts.fig = trainFig;
inferenceOpts.numBurnin = 1000;
inferenceOpts.numIter = 5000;
inferenceOpts.silent = false;
inferenceOpts.numIntensityPoints = 200;
inferenceOpts.tmin = synth.xIntensity(1);
inferenceOpts.tmax = synth.xIntensity(end);



results = evaluateInference(@runIntensityInference, synth, inferenceOpts);
results.accuracy

fig = figure(resultsFig);
clf;
fillColor = [];
yAxisMax = [];
uncertainty = true;

plotResults(results, synth, yAxisMax, fig, fillColor, uncertainty);

%% For inference over real data

% To run inference on real data (that is, non-synthetic data where there's no
% known intensity or parameters against which to evaluate inference accuracy)
% use something like

% events = getMyEventsSomehow(...);
% opts = intensityInferenceDefaults();
% % modify opts as deisred
% opts.numBurnin = 2000;
% opts.numIter = 2000;
% results = runIntensityInference(events, opts);
% plotResults(results, [], yAxisMAx, fig, fillColor, uncertainty);