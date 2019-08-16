function [x, itn] = oneProjectorMex_distNoSort(b,d,tau)
% [x, itn] = oneProjectorMex_distNoSort(b,d,tau) 
% (for distributed vectors)

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
  tau = tau/abs(d);
  [x,itn] = oneProjectorMex_I();
else
  [x,itn] = oneProjectorMex_D();
end

% Subfunctions
% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_I()
% ----------------------------------------------------------------------

    % Version without sort, currently performs best under parallelization but less tested
    % Assumes that b is a distributed array

    % Initialization
    n     = length(b);
    bNorm = gather(sum(b));
    
    % Check for quick exit.
    if (tau >= bNorm), x = b; itn = 0; return; end
    if (tau <  eps  ), x = min(0,b);  itn = 0; return; end
    
    PROJ_ITER = 0;
    relative_tol = eps(classUnderlying(b));
    
    spmd
        b_codistr = getCodistributor(b);
        local_b = getLocalPart(b);
        
        while ((bNorm - tau)/tau) > relative_tol
            alpha = (bNorm - tau) / n;
            
            %perform thresholding for each worker          
            local_b = local_b - alpha;
            nonzero_flag = local_b > 0;
            local_b = local_b .* nonzero_flag;        
            
            % compute current 1norm and number of values above 0
            local_bNorm = sum(local_b);
            local_n = sum(nonzero_flag);
            local_broadcast_vals = [local_bNorm; local_n];
            
            % broadcast to other labs
            broadcast_vals = gcat(local_broadcast_vals);
            
            % compute global 1norm and number of nonzeros
            broadcast_vals = sum(broadcast_vals,2);
            bNorm = broadcast_vals(1);
            n = broadcast_vals(2);
            
            PROJ_ITER = PROJ_ITER + 1;
        end
        
        b = codistributed.build(local_b, b_codistr, 'noCommunication');
    end
    x = b;
    
    % PROJ_ITER
    
    % Set number of iterations
    itn = gather(PROJ_ITER(1));
    
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