function y = funfun(varargin)
%% Function for testing
if nargin == 1 && varargin{1} == 0
    y = [500 300 1 1];
    
else % Multiply
    A = varargin{1};
    S = varargin{2};
    x = varargin{3};
    mode = varargin{4};
    
    if mode == 1
        y = A*x - S;
    else
        y = A'*x - conj(S);
    end
end % multiply

end % funfun