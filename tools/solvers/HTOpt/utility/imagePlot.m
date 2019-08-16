function imagePlot(z, varargin)
% IMAGEPLOT - Convenience function for creating image plots
% 
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%   imagePlot(z,'optionalarg1',optionalargval1,...)
%
% Input:
%   z - m x n real matrix to display
%
% Optional inputs
%  'cmap'          - colormap to use (default: 'default'). Will use seismic
%                    colormap if available
%  'coloraxis'     - range of the coloraxis (default: [])
%  'titleStr'      - title of the plot (default: [], no title)
%  'refVec'        - if specified, will add the relative error compared to z to the 
%                    title (default: [])
%  'cbar'          - if true, will display the colorbar (default: true)
%  'centercaxis'   - if true, will center the color axis (default: false)
%  'newFigure'     - if true, will display the results in a new figure (default: true)
%  'specPlot'      - if true, will display the 2D Fourier power spectrum of the signal z (default: false)
%  'zoomFactor'    - integer value. If > 1, will zoom in by the factor specified (default: 1)
%  'xLabel'        - label for the x axis (default: [])
%  'yLabel'        - label for the y axis (default: [])
%  'visible'       - if true, will display the resulting figure (default: true)
%                    invisible figures are useful when generating lots of images in batch mode
     
    
[dims,cmap,coloraxis,titleStr,refVec,cbar,centercaxis,newFigure,specPlot,zoomFactor,xLabel,yLabel,xAxis,yAxis,visible] = process_options(varargin, ...
                       'dims', size(z), ...
                       'cmap','default', ...
                       'coloraxis', [], ...
                       'titleStr', '', ...
                       'refVec', [], ...
                       'cbar', true, ...
                       'centercaxis',false,...
                       'newFigure', true, ...
                       'specPlot', false, ...
                       'zoomFactor', 1, ...
                       'xLabel', [], ...
                       'yLabel', [], ...
                       'xAxis', [], ...
                       'yAxis', [],...
                       'visible',true);

%Unpack option variables into workspace


if(length(dims) > 1 && dims(1)==1 && dims(2) > 1)
    dims(1) = dims(2); dims(2) = 1;
    z = z(:);
end

if(newFigure)
    if visible
        figure;
    else
        figure('Visible','Off');
    end
end

if(~isempty(refVec) && ~isreal(refVec))
    disp('imagePlot: reference vector is complex, only taking real part');
    refVec = real(refVec);
end

if(~isreal(z))
    disp('imagePlot: data is complex, only plotting real part');
    z = real(z);    
end


if(~isempty(refVec) && ~specPlot)
   relerr = norm(z(:) - refVec(:))/norm(refVec(:));
   titleStr = [titleStr, ' relerr = ' num2str(relerr)];
end

z = reshape(z,dims);
if (specPlot)
    if(size(dims,2) > 1)
        z = fftshift(abs(fft2(z)))/numel(z);
    else
        z = fftshift(abs(fft(z)))/numel(z);
    end
end
if(~isempty(xAxis) && ~isempty(yAxis))
    imagesc(xAxis, yAxis, z);
else
    imagesc(z);
end
zoom reset;
zoom(zoomFactor);

if(~isempty(coloraxis))
    caxis(coloraxis);
else
    coloraxis = caxis;
end
if(cbar)
    colorbar;
end
if(centercaxis)
    m = min(abs(coloraxis));
    caxis([-m m]);
end
if strcmp(upper(cmap),'DEFAULT')    
    if exist('seismic_colormap','file')~=0
        colormap(seismic_colormap(255));
    elseif exist('seismic','file')~=0
        colormap(seismic(255));
    else
        colormap('gray');
    end
else
    colormap(cmap);
end
title(titleStr); xlabel(xLabel); ylabel(yLabel);

end