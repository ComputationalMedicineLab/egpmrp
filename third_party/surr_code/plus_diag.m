function X = plus_diag(X, y)
%PLUS_DIAG add a scalar or vector onto the diagonal of a matrix
%
% X = plus_diag(X, y)
%
% For vector y: X = X + diag(y);
% For scalar y: X = X + diag(repmat(y, length(X), 1));
% (although more efficient code is used)
%
% Inputs:
% 	 X NxN 
% 	 y Nx1, 1xN or 1x1
%
% Outputs:
% 	 X  NxN
%
% Note: in older version of Matlab and Octave this function never adds y to X
% in place. As long as plus_diag is called from a function (not a script or the
% command-line) Matlab >= R2007a will update X in place when it can.
%
% Iain Murray, June 2006, July 2008.

[N, M] = size(X);
if N~=M, error('X must be square'), end

diagidx = (0:N-1)*N + (1:N);
X(diagidx) = X(diagidx) + y(:)';
