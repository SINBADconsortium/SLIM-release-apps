% BY : TMH 1/8 1997
% Updated by Thomas Mejer Hansen : 22-03-1999
% Updated by Tim Lin : 12-05-2010
%
function [S,nt,nx]=ReadBinFast(fileid,nt,nx,byte);
  
  if nargin==4
    if byte=='b'
      b_order='ieee-be';
    else
      b_order='ieee-le';
    end
  end
  

  if exist('b_order')==1, fid=fopen(fileid,'r',b_order);
  else, fid=fopen(fileid,'r'); end
  S=fread(fid,inf,'float');
  fclose(fid);
  
  l=length(S);
  
  if nt==0, nt=l/nx; disp(['    nt=',num2str(nt)]);end
  if nx==0, nx=l/nt; disp(['    n_trace=',num2str(nx)]); end,
  
  S=reshape(S,nt,nx);
  