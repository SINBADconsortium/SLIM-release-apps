function curvelet_denoise_gen(result_file, ratio, op)
% a generic version of curvelet_denoise without hard-coded
% parameters in the program
% input: result_file
%        ratio
%        op: struct with the following fields
%            dx, dz, nz, nx, nz_pad,
%            nbscales, nbangles, finest, is_real, ttype
%            useDepthPrec, 
%            mutestart,muteend,
dx = op.dx;
dz = op.dz;
nz = op.nz;
nx = op.nx;
nz_pad = op.nz_pad; % for dipole source

op.n = nx;
op.m = nz + nz_pad;
C = opCurvelet(op.m,op.n,op.nbscales,op.nbangles,op.finest,op.ttype,op.is_real);
C = C';
[header,nb_coeff] = get_curvelet_hdr(op);
% other operators
% real restriction op
opGetReal = opRealRestriction(op.m*op.n);
% depth preconditioner
if op.useDepthPrec
    dep = dz:dz:dz*op.m;
    opDepth = opKron(opDirac(op.n),opDiag(sqrt(dep)));
else
    opDepth = opDirac(op.n*op.m);
end
% top-muting op
opModelMute = opKron(opDirac(op.n),opLinMute(op.m,op.mutestart,op.muteend));
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