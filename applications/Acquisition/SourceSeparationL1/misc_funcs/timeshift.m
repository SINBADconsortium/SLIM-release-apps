function [ output ] = timeshift(x,params,mode)

% Time delay matrices
shiftforward_S1 = exp(-1i*params.omega*params.tdelay_S1);
shiftforward_S2 = exp(-1i*params.omega*params.tdelay_S2);
shiftadjoint_S1 = exp(1i*params.omega*params.tdelay_S1);
shiftadjoint_S2 = exp(1i*params.omega*params.tdelay_S2);

if mode == 1
   x1 = x(1:params.nr*params.nc);
   x2 = x(params.nr*params.nc+1:end);
   S1 = shiftforward_S1(:).*x1;
   S2 = shiftforward_S2(:).*x2;
   output = S1 + S2; 
else
   S1 = shiftadjoint_S1(:).*x;
   S2 = shiftadjoint_S2(:).*x;      
   output = [S1;S2];
end

end % function end

