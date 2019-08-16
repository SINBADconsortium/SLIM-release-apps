function [ varargout ] = freq_grid( varargin )
%FREQ_GRID
for i=1:length(varargin)
    xgrid = varargin{i};
    N = length(xgrid);
    dx = xgrid(2)-xgrid(1);
    fs = 1/dx;
    if mod(N,2)==0
        varargout{i} = (-N/2:(N/2-1))*fs/N;
    else
        varargout{i} = (-(N-1)/2:((N-1)/2)) *fs/N;
    end
end
end

