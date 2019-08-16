%% Matlab coder compile
%non vectorized, for-loop version
coder_compile('Helm3dmvp.m');
% vectorized version
coder_compile('Helm3dmvp_v.m');

%% Generate random inputs
n = 100*ones(3,1); npml = 10*ones(2,3); h = [25;25;25];
wn = complex(rand(prod(n),1)+0.5,rand(prod(n),1)+0.5); wn = real(wn);
q = complex(randn(prod(n),1),randn(prod(n),1));

%% Timing functions 
disp('Matrix vector product time (s)');

tic,y = Helm3dmvp(wn,h,n,npml,q,1);disp(toc);

create_helm3dmvp;
coder_compile('Helm3dmvp.m');

tic,for i=1:10,y1 = Helm3dmvp_mex(wn,h,n,npml,q,1);end,disp(toc/10);

tic,y2 = Helm3dmvp_v(wn,h,n,npml,q,1); disp(toc);

tic,y3 = Helm3dmvp_v_mex(wn,h,n,npml,q,1); disp(toc);

tic,yref = Helm3dmvp_forw_mex(wn,h,n,npml,q,8); disp(toc);
disp(norm(vec(y1)-vec(yref)));