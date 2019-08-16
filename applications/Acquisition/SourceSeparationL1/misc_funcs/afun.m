function [ output ] = afun(x, C, params, mode)

if mode == 1
   output = opBlockDiag(C',C')*vec(x);  % from curvelet domain to physical domain
   output = timeshift(output, params, 1);
else
   output = timeshift(x, params, -1);
   output = opBlockDiag(C,C)*vec(output);  % from physical domain to curvelet domain
end

end % function end

