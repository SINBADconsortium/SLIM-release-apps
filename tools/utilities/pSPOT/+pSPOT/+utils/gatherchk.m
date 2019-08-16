function y = gatherchk(tmpy,mode,gatherFlag)
% GATHERCHK     Function for checking and processing of the gathering 
%               of the result
%
%   GATHERFLAG specifies whether to gather the results to a local array
%   or leave them distributed, default is 0.
%   GATHERFLAG = 0 will leave them distributed.
%   GATHERFLAG = 1 will gather the results of forwards or adjoint multiplication.
%   GATHERFLAG = 2 will gather only in forward mode.
%   GATHERFLAG = 3 will gather only in backward (adjoint) mode.

if mode == 1
    if gatherFlag == 1 || gatherFlag == 2
        y = gather(tmpy);
    else
        y = tmpy;
    end
else % mode == 2
    if gatherFlag == 1 || gatherFlag == 3
        y = gather(tmpy);
    else
        y = tmpy;
    end
end % gather