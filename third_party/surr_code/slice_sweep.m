function particle = slice_sweep(particle, slice_fn, sigma, step_out)
%SLICE_SWEEP one set of axis-aligned slice-sampling updates of particle.pos
%
%     particle = slice_sweep(particle, slice_fn[, sigma[, step_out]])
%
% The particle position is updated with a standard univariate slice-sampler.
% Stepping out is linear (if step_out is true), but shrinkage is exponential. A
% sensible strategy is to set sigma conservatively large and turn step_out off.
% If it's hard to set a good sigma though, you should leave step_out=true.
%
% Inputs:
%     particle   sct   Structure contains:
%                              .pos - initial position on slice as Dx1 vector
%                                     (or any array)
%                           .Lpstar - log probability of .pos (up to a constant)
%                         .on_slice - needn't be set initially but is set
%                                     during slice sampling. Particle must enter
%                                     and leave this routine "on the slice".
%     slice_fn   @fn   particle = slice_fn(particle, Lpstar_min)
%                      If particle.on_slice then particle.Lpstar should be
%                      correct, otherwise its value is arbitrary.
%        sigma (D|1)x1 step size parameter(s) (default=1)
%     step_out   1x1   if non-zero, do stepping out procedure (default), else
%                      only step in (saves on fn evals, but takes smaller steps)
%
% Outputs:
%     particle   sct   particle.pos and .Lpstar are updated.

% Originally based on pseudo-code in David MacKay's text book p375
% Iain Murray, May 2004, January 2007, June 2008, January 2009

if nargin < 3; sigma = 1; end
if nargin < 4; step_out = 1; end

DD = numel(particle.pos);
if length(sigma) == 1
    sigma = repmat(sigma, DD, 1);
end

%for dd = 1:DD
% A random order is more robust generally and important inside
% algorithms like nested sampling and AIS
for dd = randperm(DD)
    Lpstar_min = particle.Lpstar + log(rand);

    % Create a horizontal interval (x_l, x_r) enclosing x_cur
    x_cur = particle.pos(dd);
    rr = rand;
    x_l = x_cur - rr*sigma(dd);
    x_r = x_cur + (1-rr)*sigma(dd);
    if step_out
        particle.pos(dd) = x_l;
        while 1
            particle = slice_fn(particle, Lpstar_min);
            if ~particle.on_slice
                break
            end
            particle.pos(dd) = particle.pos(dd) - sigma(dd);
%            fprintf('1: sigma %g particle.pos %g\n', sigma(dd), particle.pos(dd));
        end
        x_l = particle.pos(dd);
        particle.pos(dd) = x_r;
        while 1
            particle = slice_fn(particle, Lpstar_min);
            if ~particle.on_slice
                break
            end
            particle.pos(dd) = particle.pos(dd) + sigma(dd);
%            fprintf('2: sigma %g particle.pos %g (%d) \n', sigma(dd), particle.pos(dd), dd);
        end
        x_r = particle.pos(dd);
    end

    % Make proposals and shrink interval until acceptable point found
    % One should only get stuck in this loop forever on badly behaved problems,
    % which should probably be reformulated.
    chk = 0;
    while 1
        particle.pos(dd) = rand*(x_r - x_l) + x_l;
%        fprintf('x_l %g x_r %g\n', x_l, x_r);
        particle = slice_fn(particle, Lpstar_min);
        if particle.on_slice
            break % Only way to leave the while loop.
        else
            % Shrink in
            if particle.pos(dd) > x_cur
                x_r = particle.pos(dd);
            elseif particle.pos(dd) < x_cur
                x_l = particle.pos(dd);
            else
                error('BUG DETECTED: Shrunk to current position and still not acceptable.');
            end
        end
    end
end

