function op = opSaveSnapshot(n,savefilename)
% Syntax:
% op = opSaveSnapshot(n,savefilename)
%
% Description:
% Make an operator which saves each update of an inversion process.
% It is essentially an opDirac operator that saves the data. Note that the data
% is not overwritten upon update; instead, new updates are appended to the end.
%
% Input list:
% n: length of the vector you are saving.
% savefilename: file name for saving.
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

saved = zeros(n,1);

funchandle = @(x,mode) opSaveSnapshot_intrnl(n,x,mode);
op=opFunction(n,n,funchandle);

function y = opSaveSnapshot_intrnl(n,x,mode)
	switch mode  
	    case 0
	        y = {n,n,[1,1,1,1],{'opSaveSnapshot'}};
	    case 1
	        saved(:,end+1) = x(:);
            save(savefilename,'saved');
    		y = x;
	    case 2
	        y = x;
		case 3
			y = saved;
	end
end
end