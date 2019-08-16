function [x, itn] = oneProjectorMex_opt(b,d,tau)
% [x, itn] = oneProjectorMex_opt(b,d,tau) 
%
% oneProjectorMex_opt is the optimized oprtable version of oneProjectorMex 
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
%   $Id: oneProjectorMex_opt.m 37 2011-04-27 19:58:40Z timlin $
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
  tau = tau/abs(d);
  [x,itn] = oneProjectorMex_I();
else
  [x,itn] = oneProjectorMex_D();
end

% Subfunctionns
% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_I()
% ----------------------------------------------------------------------

    % Optimized portable version that eliminated for-loops

    % Initialization
    n     = length(b);
    x     = b;
    bNorm = norm(b,1);   
    
    % Check for quick exit.
    if (tau >= bNorm), x = b; itn = 0; return; end
    if (tau <  eps  ), x = min(0,b);  itn = 0; return; end
    
    b = sort(b,'descend'); % Descending.
     
    csb = cumsum(b) - tau;
    
    num = [1:n];
    csb = csb ./ num(:);
    clear num
    
    thres_index = find(csb >= b,1);
    if thres_index == 1
        alpha = 0;
    elseif isempty(thres_index)
        alpha = csb(end);
    else
        alpha = csb(thres_index-1);
    end
    
    clear b
    clear csb
    
    % Set the solution by applying soft-thresholding with
    % the previous value of alpha
    x = x - alpha;
    x = max(0,x);
    
    % Set number of iterations
    itn = thres_index;
end


% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_D()
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

end % function oneProjectorMex
