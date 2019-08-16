function op = opTukeyWinMask(total_length,alpha,mask,taper_sides)
% Syntax:
% op = opTukeyWinMask(total_length,alpha,mask,taper_both_sides)
%
% Description:
% make a 1D Tukey window. The only difference from Spot-builtin opWindow(n,
% 'Tukey', alpha) is that it admits a mask. Parameter "alpha" should be between
% 0 and 1 (including).
%
% Input list:
% total_length: length of the vector you want to apply taper to.
% alpha: ratio between the tapering-off part and the entire window. On each side
%		of the window, the ratio is alpha/2. Together it is alpha.
% mask: must be a logical vector of the same length as the window, indicating 
% 		which entries remain untouched and which are thrown away.
% taper_sides: a struct with two fields taper_sides.left and taper_sides.right, 
%		will taper the side when the field is 1.
%
% Output list:
% op: the operator
%
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Feb/14/2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if not(exist('mask','var'))
	mask = true(total_length,1);
end
mask = logical(mask);

if not(exist('taper_sides','var'))
	taper_sides.left = 1;
	taper_sides.right = 1;
end

if length(mask) ~= total_length
	error('Taper mask is not of correct size');
end

if (alpha == 0)
	taper = ones(total_length,1);
elseif (alpha == 1)
	taper = opWindow(total_length,'Hann');
	taper = diag(double(taper));
elseif ((alpha > 0) && (alpha < 1))
	taper = opWindow(total_length,'Tukey', alpha);
	taper = diag(double(taper));
else
	error('Fatal: Tapering parameter for Tukey window should be between 0 and 1.')
end

half_length = floor(total_length/2);
if (taper_sides.left == 0)
	taper(1:half_length) = taper(half_length);
end
if (taper_sides.right == 0)
	taper(half_length+1:end) = taper(half_length+1);
end

taper = full(taper(mask));

n = length(taper);

funchandle = @(x,mode) opTukeyWinMask_intrnl(x,mode);
op = opFunction(n,n,funchandle);

    function y = opTukeyWinMask_intrnl(x,mode)
        if (mode == 0)
            y = {n,n,[0,1,0,1],{'Dirac'}};
        else
            y = x.*taper(:);
        end
    end
end