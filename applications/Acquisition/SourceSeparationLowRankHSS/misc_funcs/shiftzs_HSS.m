function [ output ] = shiftzs_HSS(x, params, mode)

% Time delay matrices   
shiftforward = exp(-1i*params.omega*(params.tdelay_sub));
shiftadjoint = exp(1i*params.omega*(params.tdelay_sub));

% NOTE: params.sign is used to use entries (puts ones) inside the MH diamond, and puts zeros outside the MH diamond
if mode == 1
   x1 = x(1:params.mhnumr*params.mhnumc);
   x2 = x(params.mhnumr*params.mhnumc+1:end);
   S1 = params.sign(:).*x1;
   S2 = params.sign(:).*(shiftforward(:).*x2);
   output = S1 + S2; 
else
   S1 = params.sign(:).*x;
   S2 = params.sign(:).*(shiftadjoint(:).*x);      
   output = [S1;S2];
end

end % function end

