% ReadSuFast
%
% PURPOSE : reads a SEISMIC section i  SU format in big endian format, 
%           strips the headers and returns the field in the matrix seis.
%           If nx==0 and nt<>0, nx will be computed
%           If nt==0 and nx<>0, nt will be computed           
%
% Call : function seis=ReadSuFast(fileid,nt,nx,'byteorder');
%           byteorder : 'l' for little or 'b' for big endian (Default : Native )
%
% BY : TMH 1/8 1997
% Updated by Thomas Mejer Hansen : 22-03-1999
% Updated by Tim Lin : 12-05-2010
%
function [S,nt,nx,dt]=ReadSuFast(fileid,nt,nx,byte);
  
  if nargin==4 | nargin==2,
    if nargin==2, byte=nt; nt=0; end
    if byte=='b'
      b_order='ieee-be';
    else
      b_order='ieee-le';
    end
  else
    % USE DEFAULT BYTEORDER
  end
  
  % read first trace header to get ns and dt
  if exist('b_order')==1
      fid=fopen(fileid,'r',b_order);
  else
      fid=fopen(fileid,'r'); 
  end
  S0=fread(fid,60,'int16');
  fclose(fid);
  
  if nargin==1 | nargin==2,
    if exist('b_order')==1
        fid=fopen(fileid,'r',b_order);
    else
        fid=fopen(fileid,'r');
    end
    
    nt=S0(58);
    if nt<0,
      % THIS CANNOT BE, MAYBE WRONG BYTE ORDER
      error([mfilename,' : negative number of sample, most likeky cause is wrong endianess, try using suswapbytes'])
      S=[];nt=nt;nx=[];
      return
    end
    if nt<10,
      % THIS IS VERY UNUSUAL, MAYBE WRONG BYTE ORDER
      error([mfilename,' : very small number of sample, most likeky cause is wrong endianess, try using suswapbytes'])
      S=[];nt=nt;nx=[];
      return
    end
    disp(['Reading SU trace header, nt=',num2str(nt)])
    nx=0;
  end
  
  dt = S0(59);
  if dt<0,
    % THIS CANNOT BE, MAYBE WRONG BYTE ORDER
    error([mfilename,' : negative dt, most likeky cause is wrong endianess, try using suswapbytes'])
    S=[];nt=nt;nx=[];
    return
  end
  
  if exist('b_order')==1, fid=fopen(fileid,'r',b_order);
  else, fid=fopen(fileid,'r'); end
  S=fread(fid,inf,'float');
  fclose(fid);
  
  l=length(S);
  
  
  if nt==0, nt=l/nx-60; disp(['    nt=',num2str(nt)]);end
  if nx==0, nx=l/(nt+60); disp(['    n_trace=',num2str(nx)]); end,
  
  S=reshape(S,nt+60,nx);
  
  % remove trace headers
  S(1:60,:) = [];
  
  % return dt in seconds
  dt = dt / (1e6);
