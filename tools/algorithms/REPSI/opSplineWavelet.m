function op = opSplineWavelet(m,n, filt, smooth, levels, type)
% OPSPINEWAVELET  Wavelet operator (modified for spline wavelets)
%
%    OPSPLINEWAVELET(M,N,FILT_LENGTH,SMOOTH,LEVELS,TYPE) creates a wavelet
%    operator that acts on a 2D signal built on spline wavelets.
%   
%    M,N            - size of input domain
%    FLIT_LENGTH    - length of the filters (def 16)
%    SMOOTH         - vector with orders SMOOTH>-0.5 : degree of the fractional splines (def 2)
%                       (note this is the alpha in common fractal spline literature)
%    LEVELS         - levels used in the wavelet transform (def 5)
%    TYPE           - type of B-splines
%                     '+ortho'   (causal orthonormal)
%                     '-ortho'   (anticausal  orthonormal)  
%                     '*ortho'   (symmetric orthonormal, DEFAULT);
%                     '+bspline' (causal B-spline)
%                     '-bspline' (anticausal B-spline)
%                     '*bspline' (symmetric B-spline). 
%                     '+dual'    (causal dual)
%                     '-dual'    (anticausal dual)
%                     '*dual'    (symmetric dual).

%   Copyright 2011, Tim Lin, SLIM

if (nargin < 4), smooth = 2;   end;
if (nargin < 3), filt   = 16;   end;
if (nargin < 5), levels = 5;   end;
if (nargin < 6), type   = '*ortho';;      end;


h =  MakeFracSplinesFilters(filt,smooth,type);

subfunc_handle = @(x,mode) opWavelet_intrnl(m,n,filt,smooth,levels,type,h,x,mode);

op = opFunction(m*n,m*n,subfunc_handle); % return a SPOT operator using constructor opFunction


function y = opWavelet_intrnl(m,n,filt,smooth,levels,type,h,x,mode)
if mode == 0
   y = {n*m,n*m,[0,1,0,1],{'SplineWavelet',filt,smooth,levels,type}};
elseif mode == 2
   if isreal(x)
     [y,l] = spot.rwt.midwt(reshape(x,[m,n]), h, levels);
   else
     [y1,l] = spot.rwt.midwt(reshape(real(x),[m,n]), h, levels);
     [y2,l] = spot.rwt.midwt(reshape(imag(x),[m,n]), h, levels);
     y = y1 + sqrt(-1) * y2;    
   end
   y = reshape(y,[m*n,1]);
elseif mode == 1
   if isreal(x)
      [y,l] = spot.rwt.mdwt(reshape(x,[m,n]), h, levels);
   else
      [y1,l] = spot.rwt.mdwt(reshape(real(x),[m,n]), h, levels);
      [y2,l] = spot.rwt.mdwt(reshape(imag(x),[m,n]), h, levels);
     y = y1 + sqrt(-1) * y2;    
   end   
   y = reshape(y,[m*n,1]);
end
