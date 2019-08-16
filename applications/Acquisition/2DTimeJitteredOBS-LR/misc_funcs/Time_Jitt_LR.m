function [f1, f2] = Time_Jitt_LR(x,g,params)


% Use:
%   Time_Jitt_LR(x,g,params)
%
% Input:   
%      x - unknow data 
%      g - computed gradient       
%      params - parameter file

% Author: Rajiv Kumar
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: January, 2015

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%----------------------------------------------------------------------------------------------------


% replacement of NLfunforward (this define complete A)
Ft = opDFTR(params.nt);
Ft = opKron(opDirac(params.nr*params.nc),Ft);
L = x(1:params.nf*params.nm*params.k);
R = x(params.nf*params.nm*params.k+1:end);
L = reshape(L,params.nf,params.nm,params.k);
R = reshape(R,params.nf,params.nh,params.k);
MH = opMH(params.nr,params.nc);


if isempty(g)
    output = zeros(params.nf,params.nr*params.nc);
    for i = 1:params.nf
        output(i,:) = MH'*vec(squeeze(L(i,:,:))*squeeze(R(i,:,:))');
    end
    
    f1 = params.RM*(Ft'*vec(output));
    f2 = 0;
    clear vars output
else
    fp  = Ft*(params.RM'*vec(g));
    fp  = reshape(fp,params.nf,params.nr,params.nc);
    
    fL  = zeros(params.nf,params.nm,params.k);
    fR  = zeros(params.nf,params.nh,params.k);
    fLR = zeros(params.nf,params.nm,params.nh);
    
    for i = 1:params.nf
        ftest = reshape(MH*vec(squeeze(fp(i,:,:))),params.nm,params.nh);
        fLR(i,:,:) = ftest;
        fL(i,:,:)  = ftest*squeeze(R(i,:,:));
        fR(i,:,:)  = ftest'*squeeze(L(i,:,:));
    end
    f1 = [vec(fL); vec(fR)];
    f2 = vec(fLR);
    clear vars fp fL fR fLR
end
end
