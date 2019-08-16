function x = NormL1_project(x,weights,tau)
% x = NormL1_project(x,weights,tau)
% 
% Projects x onto a L1 ball such that norm(x,1) <= tau
% NOTE: this is a modified NormL1_project for efficient memory usage (Tim Lin, SLIM, Apr 27 2011)


% Check arguments
if nargin < 2
  error('The oneProjector function requires at least two parameters');
end
if nargin < 3
  tau = weights;
  weights   = [];
end

if isreal(x)
    x_knownToBePositive = 0;
    x_isComplex = 0;
else
    % for complex number L1 projection, we only project the absolute value of x
    x_angle = angle(x);
    x = abs(x);
    
    x_knownToBePositive = 1;
    x_isComplex = 1;
end
    
% original call is 'x  = oneProjector(x,weights,tau) .* xsign;', but oneProjector is now inlined below
% (Tim Lin, SLIM, Apr 27 2011)
% BEGIN INLINE oneProjector ============================================================
% Check weight vector
if isempty(weights), weights = 1; end;

if ~isscalar(weights) && ( length(x) ~= length(weights) )
    error('x and vector weights must have the same length');
end

% Quick return for the easy cases.
if isscalar(weights)  &&  weights == 0
    itn = 0;
    return
end

% Get sign of b and set to absolute values
if not(x_isComplex) || not(x_knownToBePositive)
    s = int8(sign(x)); % sign(x) returns a float array, but only contains {-1, 0, 1}, so cast into int
    x = abs(x);
else
    s = 1; % do nothing for s
end

% Perform the projection
if not(isa(x,'distributed'))
    if isscalar(weights)
        [x,itn] = oneProjectorMex_noSort(x,tau/weights);
    else
        weights   = abs(weights);
        idx = find(weights > eps); % Get index of all non-zero entries of d
        [x(idx),itn] = oneProjectorMex_noSort(x(idx),weights(idx),tau);
    end
else
    if isscalar(weights)
        % The NoSort version should perform best under paralleization, but accuracy might be slightly sacrificed
        [x,itn] = oneProjectorMex_distNoSort(x,tau/weights);
    else
        % d   = abs(d);
        % idx = find(d > eps); % Get index of all non-zero entries of d
        % [b(idx),itn] = oneProjectorMex_noSort(b(idx),d(idx),tau);
        error('Projection onto weighted one-norm ball for distributed vector is currently not implemented yet')
    end
end
% END INLINE oneProjector ============================================================

% Restore signs in x
x = x .* double(s);  % need to re-cast s back into floats for mult with float x
if x_isComplex
    x = x .* exp(i*x_angle);
    clear x_angle
end


end % NormL1_project()