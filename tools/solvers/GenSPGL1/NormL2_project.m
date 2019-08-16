function x = NormL2_project(x,weights,tau)
% x = NormL2_project(x,weights,tau)
% 
% Projects x onto a L2 ball such that norm(x,1) <= tau
% NOTE: this is a modified NormL2_project for efficient memory usage (Tim Lin, SLIM, Apr 27 2011)


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
    
    % original call is 'x = oneProjector(x,weights,tau);', but oneProjector is now inlined
    % (Tim Lin, SLIM, Apr 27 2011)
    oneProjector();
else
    % for complex number L1 projection, we only project the absolute value of x
        % ORIGINAL CODE:
        % xsign  = sign(x);
        % x = abs(x);
        % REPLACED WITH THE FOLLOWING TO AVOID ALLOCATING ANOTHER VECTOR THE SIZE OF X
        xamp = abs(x);
        x = sign(x);
        xsign = x;
        x = xamp;
        clear xamp
    
    x_knownToBePositive = 1;
    
    % original call is 'x  = oneProjector(x,weights,tau) .* xsign;', but oneProjector is now inlined
    % (Tim Lin, SLIM, Apr 27 2011)
    oneProjector();
    x = x .* xsign;
end


% Function oneProjector inlined to avoid creating extra copies of x
function oneProjector()
% ONEPROJECTOR  Projects x onto the weighted one-norm ball of radius tau
%
% NOTE: this function is inlined as a nested function of NormL1_project for memory efficiency
% see oneProjector.m in 'private' directory for documentation of original function
% (Tim Lin, SLIM, Apr 27 2011)

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

if(norm(x.*weights) > tau)
    x = (tau/weights).*(x/norm(x));
end

% 
% % Get sign of b and set to absolute values
% if not(x_knownToBePositive)
%     s = int8(sign(x)); % sign(x) returns a float array, but only contains {-1, 0, 1}, so cast into int
%     x = abs(x);
% else
%     s = 1; % do nothing for s
% end
% 
% % Perform the projection
% if ~isa(x,'distributed')
%     if isscalar(weights)
%       [x,itn] = oneProjectorMex_noSort(x,tau/weights);
%     else
%       weights   = abs(weights);
%       idx = find(weights > eps); % Get index of all non-zero entries of d
%       [x(idx),itn] = oneProjectorMex_noSort(x(idx),weights(idx),tau);
%     end
% else
%     if isscalar(weights)
%         % The NoSort version should perform best under paralleization, but accuracy might be slightly sacrificed
%       [x,itn] = oneProjectorMex_distNoSort(x,tau/weights);
%     else
%       % d   = abs(d);
%       % idx = find(d > eps); % Get index of all non-zero entries of d
%       % [b(idx),itn] = oneProjectorMex_noSort(b(idx),d(idx),tau);
%       error('Projection onto weighted one-norm ball for distributed vector is currently not implemented yet')
%     end
% end

% Restore signs in x
%x = x .* double(s);  % need to re-cast s back into floats for mult with float x
end % oneProjector

end % NormL1_project()