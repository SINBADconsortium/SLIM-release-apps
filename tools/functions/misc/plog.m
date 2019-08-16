% Writes the arguments both to stdout and to the log file.
% 
% Use:
%   plog(fid,arg1,arg2,arg3...)
%
% Input:
%  fid   - file identifier for the log file
%  argX  - arguments to be printed. Please, input only strings or numbers 
%          (numerical variables will be converted to string in here)
%          
% Author: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: May, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%-----------------------------------------------------------------------------
function plog(fid,varargin)

if fid==0 ; return ; end
if labindex~=1 ; return ; end

str = [];
for i=1:nargin-1
   if isa(varargin{i},'numeric')
      
      str = [str num2str(varargin{i},3)];
   else
      str = [str varargin{i}];
   end
end
fprintf(str);
if fid~=1
fprintf(fid,str);
end
end
