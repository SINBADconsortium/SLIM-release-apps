function [img,info] = GN_up(K,b,opts,img_guess)
% This function solve the linear problem Kx = b by using L1 solver
% SPGL1,
% 
% use:
% 	[img,info] = GN_up(K,b,opts,img_guess)
% 		
% Input,
% 	K: linear matrix or operator
% 	b: measurement
% 	opts: parameters for spgl1; help spgl1 in matlab for more information.
% 	img_guess: initial guess for x
% Output,
% 	img: recovered result
% 	info: a structure with optimization information
% 
% Author: Xiang Li
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: 02, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%

if strcmp(lower(opts.sptrans),'yes')
	C   = opCurvelet(opts.ny,opts.nx,opts.scale,opts.angle,opts.finst,'ME');
	A   = K*C';
else
	A  = K;
	C  = speye(size(K,2));
end


l1itr= 0;
if ~isempty(img_guess)
	x = C * img_guess(:);
else
	x = [];
end
rNorm = [];
while l1itr < opts.iterations
	[x,r,g,info]= spgl1(A,b,norm(x,1),0,x,opts);
	l1itr = l1itr + info.iter;
	rNorm = [rNorm;info.rNorm2];
end
img   = C' * x;
info.iter   = l1itr;
info.rNorm2 = rNorm;



