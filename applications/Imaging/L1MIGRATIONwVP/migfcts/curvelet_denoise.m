function curvelet_denoise(result_file, ratio)

dx = 7.62;
dz = 7.62;
nz = 500;
nx = 781;
nz_pad = 3; % for dipole source

op.n = nx;
op.m = nz + nz_pad;
op.nbscales = max(1,ceil(log2(min(op.m,op.n)) - 3));
op.nbangles = 16;
op.finest = 1;
op.is_real = 0;
op.ttype = 'ME';
C = opCurvelet(op.m,op.n,op.nbscales,op.nbangles,op.finest,op.ttype,op.is_real);
C = C';
[header,nb_coeff] = get_curvelet_hdr(op);
% other operators
% real restriction op
opGetReal = opRealRestriction(op.m*op.n);
% depth preconditioner
dep = dz:dz:dz*op.m;
opDepth = opKron(opDirac(op.n),opDiag(sqrt(dep)));
% top-muting op
opModelMute = opKron(opDirac(op.n),opLinMute(op.m,90,100));
% concatenate operators
opRight = opModelMute*opDepth*opGetReal*C;

load(result_file,'x')

[temp, temp_idx] = sort(abs(x),'descend');
ratio_temp = sqrt(cumsum(temp.^2))/norm(temp);
temp_idx(find(ratio_temp > (1-ratio))) = [];
x_res = zeros(length(x),1);
x_res(temp_idx) = x(temp_idx);
dm = opRight*x_res;
dm = reshape(dm, op.m, op.n);
dm = dm(nz_pad+1:end,:);
noise = opRight*(x-x_res);
noise = reshape(noise, op.m, op.n);
noise = noise(nz_pad+1:end,:);

new_result_file = [result_file(1:end-4),'_denoised.mat'];
save(new_result_file,'dm','noise')