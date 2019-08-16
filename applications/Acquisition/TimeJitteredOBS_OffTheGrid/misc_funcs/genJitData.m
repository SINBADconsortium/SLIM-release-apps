function newD = genJitData(nboats, jitacq, D, ds)

%------------------------------------------------------------------------------
% genJitData generates a data cube with jittered shots
% 
% Use:
%   genJitData(nboats, jitacq, D, ds)
%
% Input: 
%     nboats - number of boats (or source vessels)
%     jitacq - a structure array including the jittered acquisition parameters 
%              (output of the function: jitter_airgunarraysJune)
%          D - input data cube (regularly sampled)     
%         ds - underlying (regular) source sampling interval

% Author: Haneet Wason
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmospheric Sciences
%         The University of British Columbia
%         
% Date: June, 2013

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%------------------------------------------------------------------------------


% Size of input data
nt = size(D,1);
nr = size(D,2);
ns = size(D,3);

% Regular trace locations for every receiver gather
regloc = linspace(0, 1-1/ns, ns)' - 0.5;

% Jittered trace locations
if nboats == 1
   jitpos = sort([jitacq.sjitb1arr1 jitacq.sjitb1arr2]);
elseif nboats == 2
   jitpos = sort([jitacq.sjitb1arr1 jitacq.sjitb1arr2 jitacq.sjitb2arr1 jitacq.sjitb2arr2]);
end
jitloc = (jitpos/(ns*ds)) - 0.5;

% Fourier operators
NFreg = opNFFT(ns, regloc);
NF = opNFFT(ns, jitloc);

% Fourier transform using regular trace locations
x = zeros(nt,nr,ns);
for k = 1:nr
    data = squeeze(D(:,k,:));
    for i = 1:nt 
        x(i,k,:) = NFreg'*(data(i,:).'); 
    end
end

% Generate jittered shots
newD = zeros(nt,nr,length(jitloc));
for k = 1:nr
    rec = squeeze(x(:,k,:));
    for i = 1:nt 
        newD(i,k,:) = NF*(rec(i,:).');     
    end
end

clear data rec

end  % function end

