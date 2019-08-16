function [out,residual_y]=apply_sys_mat(in,WtW_inv,W,analysis)

out=(W')*in;
out=WtW_inv*out;
out=W*out;
out=in-out;

if analysis==1
    residual_y=W'*out;
    residual_y=W*residual_y;
    residual_y=(residual_y+out)-in;
else
    residual_y=inf;
end

end
%sys_mat_old = (speye(size(V,1)) + HtPt*HtPt')*Pinv_old;