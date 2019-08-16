function curvelet_denoise_wCurvSyn(result_file, ratio)

dx = 7.62;
dz = 7.62;
nz = 500;
nx = 781;

op.n = nx;
op.m = nz;
op.nbscales = max(1,ceil(log2(min(op.m,op.n)) - 3));
op.nbangles = 16;
op.finest = 1;
op.is_real = 0;
op.ttype = 'ME';
C = opCurvelet(op.m,op.n,op.nbscales,op.nbangles,op.finest,op.ttype,op.is_real);
C = C';

new_result_file = [result_file(1:end-4),'_denoised.mat'];
if not(exist(new_result_file,'file'))
    load(result_file,'dm')
    x = spgl1(C, dm(:), 0, 1e-3*norm(dm(:),2));
    save(new_result_file, 'x')
else
    load(new_result_file, 'x')
end

[temp, temp_idx] = sort(abs(x),'descend');
ratio_temp = sqrt(cumsum(temp.^2))/norm(temp);
temp_idx(find(ratio_temp > (1-ratio))) = [];
x_res = zeros(length(x),1);
x_res(temp_idx) = x(temp_idx);
dm = C*x_res;
dm = reshape(dm, op.m, op.n);
noise = C*(x-x_res);
noise = reshape(noise, op.m, op.n);

save(new_result_file,'dm','noise','-append')