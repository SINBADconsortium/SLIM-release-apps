function op = opTrimBorder(nz, nx, nbz, nbx)
% function op = opTrimBorder(nz, nx)
% Trim border of the image in each iteration

% create function handle
subfunc_handle = @(x,mode) opTrimBorder_intrnl(x,mode);
m = nz*nx;
n = m;

op = opFunction(m,n,subfunc_handle); % return a SPOT operator using constructor opFunction


    function y = opTrimBorder_intrnl(x,mode)
        if mode == 0
            y = {m,n,[0,1,0,1],{'opTrimBorder'}};
        else
            y = reshape(x,nz,nx);
            y(:,1) = y(:,1)/sqrt(nbx+1);
            y(:,nx) = y(:,nx)/sqrt(nbx+1);
            y(nz,:) = y(nz,:)/sqrt(nbz+1);
            y = y(:);
        end
    end
end