function y = funBlockDiag(varargin)
%% Function for testing
if nargin == 1 && varargin{1} == 0
    y = [500 300 1 1];
    
else % Multiply
    A = varargin{1};
    x = varargin{2};
    
    y = A*x;
end % multiply

end % funBlockDiag