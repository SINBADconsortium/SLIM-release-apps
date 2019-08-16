function x = matldivide(op,b,mode)
%\  Matlab built-in left matrix divide

%   X = bimldivide(op,b) is Matlab's builtin backslash operator,
%   except that op is always a Spot operator, and b is always a numeric column
%   vector.

% if mode == 1
%     x = builtin('mldivide',double(op),b);
% else
%     x = builtin('mldivide',double(op'),b);
% end

x = lsqrdivide(op,b,mode);