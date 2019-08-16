function y = curvelet_angle_res(x,m,n,s,a,sx,sz)
% x: image.
% m,n: size of x;
% s,a: scale and angle of curvelet transform.
% sp times: default 2.
if nargin< 7, sz = 0;end
if nargin< 6, sx = 2;end

isrl = isreal(x);

 x = reshape(x,m,n);
 cc = mefcv2(x,m,n,s,a);
 aa = cc{1};
 bb = aa{1};
 [mm,nn] = size(bb);
 P = opKron(opSmooth(nn,sx),opSmooth(mm,sz));
 aa{1} = reshape(P * bb(:),size(bb));
 cc{1} = aa;
 
 y = meicv2(cc,m,n,s,a);

 if isrl
     y = real(y);
 end
 
 
 
 
 