function op = opSplineWaveletSPOT(m, n, filt, smooth, levels, type)

%------------------------------------------------------------------------------
% OPSPINEWAVELETSPOT Wavelet operator (modified for spline wavelets)
%
%   OPWAVELET(M,N,FAMILY,FILTER,LEVELS,TYPE) creates a wavelet
%   operator of given FAMILY, for M by N matrices. The wavelet
%   transformation is computed using the Rice Wavelet Toolbox.
%
%   The remaining three parameters are optional. FILTER = 8
%   specifies the filter length and must be even. LEVELS = 5 gives
%   the number of levels in the transformation. Both M and N must
%   be divisible by 2^LEVELS. TYPE = 'min' indictates what type of
%   solution is desired; 'min' for minimum phase, 'max' for
%   maximum phase, and 'mid' for mid-phase solutions. 

%   Copyright 2007, Rayan Saab, Ewout van den Berg and Michael P. Friedlander
%   http://www.cs.ubc.ca/labs/scl/sparco
%   $Id: opWavelet.m 577 2007-09-17 01:48:28Z mpf $

%   April 2013 - Wrapped into a SPOT operator
%------------------------------------------------------------------------------

if (nargin < 3), filt   = 128; end;
if (nargin < 4), smooth = 3; end;
if (nargin < 5), levels = 5; end;
if (nargin < 6), type   = '*ortho'; end;

h = MakeFracSplinesFilters(filt, smooth, type);
fh = @(x,mode) opWavelet_intrnl(m, n, filt, smooth, levels, type, h, x, mode);
op = opFunction(m*n, n*m, fh);

%==============================================================================%

function y = opWavelet_intrnl(m, n, filt, smooth, levels, type, h, x, mode)

if mode == 0
   y = {n*m, n*m, [0,1,0,1], {'SplineWaveletSPOT', filt, smooth, levels, type}};   
elseif mode == 1
   if isreal(x)
     [y,l] = spot.rwt.midwt(reshape(x,[m,n]), h, levels);
   else
     [y1,l] = spot.rwt.midwt(reshape(real(x),[m,n]), h, levels);
     [y2,l] = spot.rwt.midwt(reshape(imag(x),[m,n]), h, levels);
     y = y1 + sqrt(-1) * y2;    
   end
   y = reshape(y,[m*n,1]);   
else
   if isreal(x)
      [y,l] = spot.rwt.mdwt(reshape(x,[m,n]), h, levels);
   else
      [y1,l] = spot.rwt.mdwt(reshape(real(x),[m,n]), h, levels);
      [y2,l] = spot.rwt.mdwt(reshape(imag(x),[m,n]), h, levels);
     y = y1 + sqrt(-1) * y2;    
   end   
   y = reshape(y,[m*n,1]);   
end

end  % internal function end

end  % external function end

