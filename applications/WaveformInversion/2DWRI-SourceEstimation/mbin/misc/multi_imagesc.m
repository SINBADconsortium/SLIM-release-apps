function multi_imagesc(opts,D,varargin)
% MULTI_IMAGESC - Produces multiple identically scaled/labelled imagesc plots, perfect for all your interpolation results!
% 
%  Curt Da Silva, 2015
%  curtd@math.ubc.ca
%
%  Usage:
%     multi_imagesc(opts,D,{D1,D2,...});
%  
%  Input:
%     opts         - struct containing plotting options
%     available options:
%       .xlabel, .ylabel - x, y axis names
%       .caxis           - color axis to use (default: caxis corresponding to D)
%       .caxis_ctr       - if true, center color axis
%       .axes_fontsize   - font size of axis labels
%       .axes_fontname   - font of axis labels
%       .axes_fontweight - 'bold',etc.
%       .cmap            - colormap to use ('seismic' or 'gray' is standard)
%       .xtick, .ytick   - location of label ticks to use for the x-, y-axes
%       .xticklabel,     - name of label ticks to use for the x-, y-axes
%         .yticklabel 
%
%     D            - reference matrix to plot
%     {D1,D2, ...} - optional matrices that are plotted according to the same specifications as D
% 
%  Example:
%    btrue - true matrix, bsub - subsampled matrix, best - estimated matrix
%    opts = struct;
%    opts.xlabel = 'receiver x'; opts.ylabel = 'receiver y';
%    opts.caxis_ctr = true; 
%    opts.axes_fontsize = 24;
%    opts.axes_fontname = 'helvetica';
%    opts.axes_fontweight = 'bold';
%    opts.cmap = cmap;
%    opts.ytick = [1:10:50];
%    opts.yticklabel = {'1';'10';'20';'30';'40';'50'};
%    opts.xtick = [1:10:50];
%    opts.xticklabel = {'1';'10';'20';'30';'40';'50'};
%    close all;
%    multi_imagesc(opts,btrue,bsub,best,btrue-best);
%    saveNames = {'true','sub',['est-snr-' num2str(SNR(vec(b),vec(best)),3)],'diff'};
%    for i=1:4
%      figure(i); print(gcf,'-depsc',[saveFigsDir saveNames{i} '.eps']);
%    end

isOptionTrue = @(x) (isfield(opts,x)==1) && (opts.(x)==1);

if isfield(opts,'subplot'),subplots = opts.subplot; else subplots = false; end
nplots = 1+length(varargin);    
horig = figure; 
if subplots, subplot(1,nplots,1); end
imagesc(D); 
axorig = gca; trueaxis = caxis;
if isfield(opts,'caxis'), caxis(opts.caxis); trueaxis = caxis; end
if isOptionTrue('caxis_ctr'), trueaxis(2) = -trueaxis(1); caxis(trueaxis); end
if isfield(opts,'cmap'), colormap(opts.cmap); else colormap('jet'); end


axOpts = {}; 
in_axes_params = {'axes_fontsize','axes_fontname','axes_fontweight','ytick','xtick','yticklabel','xticklabel'};
out_axes_params = {'FontSize','FontName','FontWeight','YTick','XTick','YTickLabel','XTickLabel'};

if isfield(opts,'xlabel'), xlabel(opts.xlabel); end
if isfield(opts,'ylabel'), ylabel(opts.ylabel); end
for i=1:length(in_axes_params)
    if isfield(opts,in_axes_params{i}), 
        set(gca,out_axes_params{i},opts.(in_axes_params{i}));
        axOpts{end+1} = out_axes_params{i}; axOpts{end+1} = opts.(in_axes_params{i});       
    end   
end

% Display each figure with the same properties as the first figure
for i=1:length(varargin)
    if subplots
        subplot(1,nplots,i+1);
    else               
        hnew = figure; 
    end
    imagesc(varargin{i});
    if isfield(opts,'xlabel'), xlabel(opts.xlabel); end
    if isfield(opts,'ylabel'), ylabel(opts.ylabel); end    
    for j=1:length(out_axes_params)
        set(gca,out_axes_params{j},get(axorig,out_axes_params{j}));
    end
    if isfield(opts,'cmap'), colormap(opts.cmap); else colormap('jet'); end
    caxis(trueaxis);
end
% If using multiple subplots, make sure the figure is properly sized
if subplots,
   pos = get(horig,'Position');
   pos(3) = 1.2*nplots*pos(4);
   set(horig,'Position',pos);
end



end