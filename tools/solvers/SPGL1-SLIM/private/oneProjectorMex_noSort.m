function [x, itn] = oneProjectorMex_noSort(x,d,tau)
% [x, itn] = oneProjectorMex_noSort(b,d,tau) 
%
% oneProjectorMex_noSort is the a version of oneProjectorMex that does not rely on sorting
% (The main feature is that this method is much more readily parallelizeable)
%
% Return the orthogonal projection of the vector b >=0 onto the
% (weighted) L1 ball. In case vector d is specified, matrix D is
% defined as diag(d), otherwise the identity matrix is used.
%
% On exit,
% x      solves   minimize  ||b-x||_2  st  ||Dx||_1 <= tau.
% itn    is the number of elements of b that were thresholded.
%
% See also spgl1, oneProjector, oneProjectorMex.


%   oneProjectorMex_opt.m
%   $Id: oneProjectorMex_opt.m 20 2010-08-09 16:41:22Z timlin $
%
%   ----------------------------------------------------------------------
%   This file is part of SPGL1 (Spectral Projected Gradient for L1).
%
%   Copyright (C) 2007 Ewout van den Berg and Michael P. Friedlander,
%   Department of Computer Science, University of British Columbia, Canada.
%   All rights reserved. E-mail: <{ewout78,mpf}@cs.ubc.ca>.
%
%   SPGL1 is free software; you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as
%   published by the Free Software Foundation; either version 2.1 of the
%   License, or (at your option) any later version.
%
%   SPGL1 is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
%   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General
%   Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with SPGL1; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
%   USA
%   ----------------------------------------------------------------------

if nargin < 3
   tau = d;
   d   = 1;
end

if isscalar(d)
  [x,itn] = oneProjectorMex_I(x,tau/abs(d));
else
  [x,itn] = oneProjectorMex_D(x,d,tau);
end

end % function oneProjectorMex


% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_I(x,tau)
% ----------------------------------------------------------------------

    % Version without sort, currently performs best under parallelization but less tested

    % Initialization
    n     = length(x);
    bNorm = sum(x);   
    
    % Check for quick exit.
    if (tau >= bNorm), itn = 0; return; end
    if (tau <  eps  ), x = min(0,x);  itn = 0; return; end
    
    PROJ_ITER = 0;
    relative_tol = eps(class(x));
    
    while (real(bNorm - tau)/tau) > relative_tol
         
        alpha = (bNorm - tau) / n;
        % Set the solution by applying soft-thresholding with
        % the previous value of alpha
        x = x - alpha;
        nonzero_flag = x > 0;
        x = x .* nonzero_flag;
        n = sum(nonzero_flag);
        clear nonzero_flag
        
        bNorm = sum(x); 
        PROJ_ITER = PROJ_ITER + 1;
    end    
    % PROJ_ITER
    
    % Set number of iterations
    itn = PROJ_ITER;
end


% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_D(b,d,tau)
% ----------------------------------------------------------------------

   % Initialization
   n = length(b);
   x = zeros(n,1);

   % Check for quick exit.
   if (tau >= norm(d.*b,1)), x = b; itn = 0; return; end
   if (tau <  eps         ),        itn = 0; return; end

   % Preprocessing (b is assumed to be >= 0)
   [bd,idx] = sort(b ./ d,'descend'); % Descending.
   b  = b(idx);
   d  = d(idx);

   % Optimize
   csdb = 0; csd2 = 0;
   soft = 0; alpha1 = 0; i = 1;
   while (i <= n)
      csdb = csdb + d(i).*b(i);
      csd2 = csd2 + d(i).*d(i);
  
      alpha1 = (csdb - tau) / csd2;
      alpha2 = bd(i);

      if alpha1 >= alpha2
         break;
      end
    
      soft = alpha1;  i = i + 1;
   end
   x(idx(1:i-1)) = b(1:i-1) - d(1:i-1) * max(0,soft);

   % Set number of iterations
   itn = i;
end
