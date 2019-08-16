function [hdr,cn] = get_curvelet_hdr(op)
if strcmp(op.ttype,'ME')
    C = mefcv2(randn(op.m,op.n),op.m,op.n,op.nbscales,op.nbangles);

    hdr{1}{1} = size(C{1}{1});
    cn = prod(hdr{1}{1});
    for i = 2:op.nbscales
        nw = length(C{i});
        hdr{i}{1} = size(C{i}{1});
        hdr{i}{2} = size(C{i}{nw/2+1});
        cn = cn + nw/2*prod(hdr{i}{1}) + nw/2*prod(hdr{i}{2});
    end
else
    [tmphdr, cn] = fdct_sizes_mex(op.m,op.n,op.nbscales,op.nbangles,logical(op.finest));
    hdr = cell(1,op.nbscales);
    hdr{1} = {[tmphdr{1:2}]}; 
    for i = 2:op.nbscales - (~op.finest)
        j = 3 + 5*(i-2);
        hdr{i}={[tmphdr{j+1:j+2}];[tmphdr{j+3:j+4}];[tmphdr{j}]};
    end
    if ~op.finest,  
        hdr{end} = {[tmphdr{end-1:end}];1}; 
    end;
end