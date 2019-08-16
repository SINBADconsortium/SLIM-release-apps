function sol = analytical_wavefield(n,c_vel,h,f,s)

vel = ones(n)*c_vel;

% Analytical solution of the homogeneous wave equation
%-----------------------------------------------------
omega = 2*pi*f;
sol = zeros(n(1),n(2),n(3));
cmplx = sqrt(-1);
for k=1:n(3)
   for j=1:n(2)
      for i=1:n(1)
         x = norm(s - ([i  j  k]-1)*h);
         wn = omega/vel(i,j,k);        % For velocity
%              wn = omega*sqrt(vel(i,j,k));  % For slowness square
         sol(i,j,k) = exp(cmplx*wn*x)/(4*pi*x);
      end
   end
end

sol = real(sol);
end
