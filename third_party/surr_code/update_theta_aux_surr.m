function [theta, ff, aux, cholK] = update_theta_aux_surr(theta, ff, Lfn, Kfn, aux, theta_Lprior, slice_width)
%UPDATE_THETA_AUX_SURR MCMC update to GP hyper-param based on aux. noisy vars
%
%     [theta, ff] = update_theta_aux_noise(theta, ff, Lfn, Kfn, aux, theta_Lprior);
%
% Inputs:
%             theta Kx1 hyper-parameters (can be an array of any size)
%                ff Nx1 apriori Gaussian values
%               Lfn @fn Log-likelihood function, Lfn(ff) returns a scalar
%               Kfn @fn Kfn(theta) returns NxN covariance matrix
%                       NB: this should contain jitter (if necessary) to
%                       ensure the result is positive definite.
%
% Specify aux in one of three ways:
% ---------------------------------
%           aux_std Nx1 std-dev of auxiliary noise to add to each value
%                       (can also be a 1x1).
% OR
%        aux_std_fn @fn Function that returns auxiliary noise level(s) to use:
%                       aux_std = aux_std_fn(theta, K);
% OR
%               aux cel A pair: {aux_std_fn, aux_cache} called like this:
%                       [aux_std, aux_cache] = aux_std_fn(theta, K, aux_cache);
%                       The cache could be used (for example) to notice that
%                       relevant parts of theta or K haven't changed, and
%                       immediately returning the old aux_std.
% ---------------------------------
%
%      theta_Lprior @fn Log-prior, theta_Lprior(theta) returns a scalar
%
% Outputs:
%             theta Kx1 updated hyper-parameters (Kx1 or same size as inputted)
%                ff Nx1 updated apriori Gaussian values
%               aux  -  Last aux_std computed, or {aux_std_fn, aux_cache},
%                       depending on what was passed in.
%             cholK NxN chol(Kfn(theta))
%
% The model is draw g ~ N(0, K + S), (imagine as f ~ N(0, K) + noise with cov S)
% Draw f ~ N(m_p, C_p), using posterior mean and covariance given g.
% But implement that using nu ~ randn(N,1). Then clamp nu's while changing K.
%
% K is obtained from Kfn.
% S = diag(aux_std.^2), or for scalar aux_std (aux_std^2 * eye(N)).

% Iain Murray, November 2009, January 2010, May 2010

% If there is a good reason for it, there's no real reason full-covariance
% auxiliary noise couldn't be added. It would just be more expensive as sampling
% would require decomposing the noise covariance matrix. For now this code
% hasn't implemented that option.

% Slight modification by Tom Lasko, September 2013, to attempt recovery from
% poorly chosen auxiliary noise, as noted below (marked with TL).

N = numel(ff);

% Start constructing the struct that will be passed around while slicing
pp = struct('pos', theta, 'Kfn', Kfn);
if isnumeric(aux)
    % Fixed auxiliary noise level
    pp.adapt_aux = 0;
    pp.aux_std = aux;
    pp.aux_var = aux.*aux;
elseif iscell(aux)
    % Adapting noise level, with computations cached
    pp.adapt_aux = 2;
    pp.aux_fn = aux{1};
    pp.aux_cache = aux{2};
else
    % Simple function to choose noise level
    pp.adapt_aux = 1;
    pp.aux_fn = aux;
end
pp.gg = zeros(N, 1);
pp = theta_changed(pp);

% Instantiate g|f
pp.gg = ff(:) + randn(N, 1).*pp.aux_std;
pp.Sinv_g = pp.gg ./ pp.aux_var;

% Instantiate nu|f,gg
pp.nu = pp.U_invR*ff(:) - pp.U_invR'\pp.Sinv_g;

% Compute current log-prob (up to constant) needed by slice sampling:
theta_unchanged = true; % theta hasn't moved yet, don't recompute covariances
pp = eval_particle(pp, -Inf, Lfn, theta_Lprior, theta_unchanged);

% Slice sample update of theta|g,nu
step_out = (slice_width > 0);
slice_width = abs(slice_width);
slice_fn = @(pp, Lpstar_min) eval_particle(pp, Lpstar_min, Lfn, theta_Lprior);
pp = slice_sweep(pp, slice_fn, slice_width, step_out);
theta = pp.pos;
% debugging only:
%theta
ff = reshape(pp.ff, size(ff));

if iscell(aux)
    aux = {pp.aux_fn, pp.aux_cache};
else
    aux = pp.aux_std;
end
cholK = pp.U;
end

function pp = theta_changed(pp)
% Will call after changing hyperparameters to update covariances and
% their decompositions.
theta = pp.pos;
K = pp.Kfn(theta);
if pp.adapt_aux
    if pp.adapt_aux == 1
        pp.aux_std = pp.aux_fn(theta, K);
    elseif pp.adapt_aux == 2
        [pp.aux_std, pp.aux_cache] = pp.aux_fn(theta, K, pp.aux_cache);
    end
    pp.aux_var = pp.aux_std .* pp.aux_std;
    pp.Sinv_g = pp.gg ./ pp.aux_var;
end
pp.U = chol(K);

pp.iK = inv(K);
try % TL: try/catch to attempt recovery if aux_var is inappropriate. This is 
    % a total hack - it would be more correct to use the 
    % aux_std_fn machinery provided for in the original design. 
    pp.U_invR = chol(plus_diag(pp.iK, 1./pp.aux_var));
catch err
    fprintf('Error decomposing matrix:\n%s\nTrying to recover using larger added diagonal...', err.message)
    fprintf('\nCond %f\n',cond(plus_diag(pp.iK, 1./pp.aux_var)));
    
    % TL: Make one attempt at recovery by adding larger diagonal, then crash if
    % insufficient.
    eps = 1e-5 * max(diag(pp.iK));
    fprintf('Trying matrix with cond %f...\n', cond(plus_diag(pp.iK, eps)));
    pp.U_invR = chol(plus_diag(pp.iK, eps));
    fprintf('Done.\n');
end
%pp.U_noise = chol(plus_diag(K, aux_var_vec));
end

function pp = eval_particle(pp, Lpstar_min, Lfn, theta_Lprior, theta_unchanged)

% Prior on theta
theta = pp.pos;
Ltprior = theta_Lprior(theta);
if Ltprior == -Inf
    pp.on_slice = false;
    return;
end
if ~exist('theta_unchanged', 'var') || (~theta_unchanged)
    pp = theta_changed(pp);
end

% Update f|gg,nu,theta
pp.ff = pp.U_invR\pp.nu + solve_chol(pp.U_invR, pp.Sinv_g);

% Compute joint probability and slice acceptability.
% I have dropped the constant: -0.5*length(pp.gg)*log(2*pi)
%Lgprior = -0.5*(pp.gg'*solve_chol(pp.U_noise, pp.gg)) - sum(log(diag(pp.U_noise)));
%pp.Lpstar = Ltprior + Lgprior + Lfn(pp.ff);
%
% This version doesn't need U_noise, but commenting out the U_noise line and
% using this version doesn't actually seem to be faster?
Lfprior = -0.5*(pp.ff'*solve_chol(pp.U, pp.ff)) - sum(log(diag(pp.U)));
%fprintf('Lfprior %g\n', Lfprior);
LJacobian = -sum(log(diag(pp.U_invR)));
%fprintf('LJacobian %g\n', LJacobian);
%LJacobian = sum(log(diag(pp.U_R)));
Lg_f = -0.5*sum((pp.gg - pp.ff).^2)./pp.aux_var - sum(log(pp.aux_var.*ones(size(pp.ff))));
%fprintf('Lg_f %g\n', Lg_f);
%fprintf('Lfn(pp.ff) %g\n', Lfn(pp.ff));
pp.Lpstar = Ltprior + Lg_f + Lfprior + Lfn(pp.ff) + LJacobian;
%fprintf('pp.Lpstar %g min %g\n', pp.Lpstar, Lpstar_min);
pp.on_slice = (pp.Lpstar >= Lpstar_min);
end

