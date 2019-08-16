function [x, itn] = oneProjectorMex_distIter(b,d,tau)
% [x, itn] = oneProjectorMex_distIter(b,d,tau) 
% Return the orthogonal projection of the vector b >=0 onto the
% (weighted) L1 ball. In case vector d is specified, matrix D is
% defined as diag(d), otherwise the identity matrix is used.
%
% This version is modified to work on distributed vectors
% introduced in Matlab PCT v4.2 (2009b)
% 
% On exit,
% x      solves   minimize  ||b-x||_2  st  ||Dx||_1 <= tau.
% itn    is the number of elements of b that were thresholded.
%
% See also spgl1, oneProjector.

%   oneProjectorMex.m
%   $Id: oneProjectorMex_distIter.m 7 2010-05-27 03:57:59Z timlin $
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
  [x,itn] = oneProjectorMex_I(b,tau/abs(d));
else
  [x,itn] = oneProjectorMex_D(b,d,tau);
end

end % function oneProjectorMex


% ----------------------------------------------------------------------
function [x,itn] = oneProjectorMex_I(b,tau)
% ----------------------------------------------------------------------

   % Initialization
   n     = length(b);
   x     = b;
   bNorm = norm(b,1);
   inCore_ratio = 0.1; % vectors of this fraction of original x can be safely brought local
   
   
   % Check for quick exit.
   if (tau >= bNorm), x = b; itn = 0; return; end
   if (tau <  eps  ), x = min(0,b);  itn = 0; return; end
   
   % Modified for MATLAB PCT: iterative divide-n-conquer algorithm to sort partial array b
   
   split_b = cell(15,1);  % split_b will act as a stack, the "top" corresponds to the largest value
   split_point = cell(15,1);

   split_b{1} = b;

   ratio = 1;
   k = 1;
   while ratio > inCore_ratio

       	prev_split_b = split_b{k};
       	split_point{k} = mean(prev_split_b);
       	split_b{k+1} = prev_split_b( prev_split_b >= split_point{k} );

       	ratio = length(split_b{k+1}) / n;

       	k = k+1;
   end
   
   count = 1;
   % gotAlpha = 0;
   csb_prev      = -tau;
   % alphaPrev = 0;
   
   while 1 % keep searching inside current vector for threshold value, otherwise go one level out
       % Preprocessing (b is assumed to be >= 0)
       
       b_local = undist(split_b{k});
       b_local = sort(b_local,'descend'); % Descending.
       
       csb = cumsum(b_local) + csb_prev;
       csb_last = csb(end);
       
       num = [count:count+length(b_local)-1];
       csb = csb ./ num(:);
       clear num

       thres_index = find(csb>b_local,1) - 1;
       
       if ~isempty(thres_index)
           if thres_index < 1
               if count == 1
                   alpha = 0;
               else
                   alpha = csb_prev / (count - 1);
               end
           else
               alpha = csb(thres_index);
           end
           count = count + thres_index;
           break;
       end
       
       csb_prev = csb_last;
       count = count + length(b_local);
       
       k = k - 1;
       if (k == 0) 
           disp('WARNING: looked inside the entire vector');
           alpha = csb_prev / n;
           break;
       end
       
       next_b = split_b{k};
       split_b{k} = next_b(next_b < split_point{k});
       
       while (length(split_b{k}) / n) > inCore_ratio
           % safeguard for bringing too_large vectors into memory
           
           prev_split_b = split_b{k};
           split_point{k} = mean(prev_split_b);
           if split_point{k} == 0 % guard against all vector zeros
               % remaining vectors are zero, short circuit
               k = 1;
               split_b{k} = zeros(10,1);
               break;
           end
           split_b{k+1} = prev_split_b(prev_split_b > split_point{k});
           k = k+1;
           % disp('in bottom portion')
       end
   end
       
   % Set the solution by applying soft-thresholding with
   % the previous value of alpha
   x = max(0,x - alpha);

   % Set number of iterations
   itn = count;
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
