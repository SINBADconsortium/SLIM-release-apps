function vel_plot(v,model,opts)
% Plot velocity model
% vel_plot(v,model,opts)
%
% opts.
%    cmap      - colormap (default: 'default')
%    cbar      - show colorbar (default: false)
%    title_str - title string (default: [])
%    cax       - reference color axis (default: [])
%    slice_dir - for 3D models, direction to slice, one of 'x','y','z'
%    slice_depth - for 3D models, depth of slice
%
    v = reshape(v,model.n);
    if isfield(opts,'cmap'), cmap = opts.cmap; else cmap = 'default'; end
    if isfield(opts,'cbar'), cbar = opts.cbar; else cbar = false; end
    if isfield(opts,'title_str'),title_str = opts.title_str; else title_str = ''; end
    if isfield(opts,'cax'), cax = opts.cax; else cax = []; end
    if isfield(opts,'new_fig'), new_fig = opts.new_fig; else new_fig = false; end
    if length(model.n)==3 && model.n(3)>1
        mode3d = true;
    else
        mode3d = false;
    end
    if mode3d
        assert(isfield(opts,'slice_dir'),'Need slice_dir');
        
        if isfield(opts,'slice_idx')
            slice_idx = opts.slice_idx;
        else
            slice_idx = [];
        end
        zax = round(model.d(3)/1e3 *(0:model.n(3)-1),1);
        yax = round(model.d(2)/1e3 *(0:model.n(2)-1),1);
        xax = round(model.d(1)/1e3 *(0:model.n(1)-1),1);
        npts = 5;
        
        switch opts.slice_dir
          case 'x'
            if isempty(slice_idx)
                slice_idx = floor(opts.slice_depth/model.d(1));
            end
            slice = @(x) fliplr(rot90(squeeze(x(slice_idx,:,:)),-1));
            xlbl = 'y [km]'; ylbl = 'z [km]';
            
            xtick = floor(linspace(1,model.n(2),npts));            
            xticklabel = cellfun(@num2str,num2cell(yax(xtick)),'UniformOutput',false);
            ytick = floor(linspace(1,model.n(3),npts));            
            yticklabel = cellfun(@num2str,num2cell(zax(ytick)),'UniformOutput',false);
            
          case 'y'
            if isempty(slice_idx)
                slice_idx = floor(opts.slice_depth/model.d(2));
            end
            slice = @(x) rot90(squeeze(x(:,slice_idx,:)),-1);
            xlbl = 'x [km]'; ylbl = 'z [km]';
            xtick = floor(linspace(1,model.n(1),npts));            
            xticklabel = cellfun(@num2str,num2cell(xax(xtick)),'UniformOutput',false);
            ytick = floor(linspace(1,model.n(3),npts));            
            yticklabel = cellfun(@num2str,num2cell(zax(ytick)),'UniformOutput',false);

          case 'z'
            if isempty(slice_idx)
                slice_idx = floor(opts.slice_depth/model.d(3));
            end
            slice = @(x) squeeze(x(:,:,slice_idx));
            xlbl = 'x [km]'; ylbl = 'y [km]';
            xtick = floor(linspace(1,model.n(1),npts));            
            xticklabel = cellfun(@num2str,num2cell(xax(xtick)),'UniformOutput',false);
            ytick = floor(linspace(1,model.n(2),npts));            
            yticklabel = cellfun(@num2str,num2cell(yax(ytick)),'UniformOutput',false);

          otherwise
            error('Invalid slice_dir');
        end
    end
    if new_fig
        figure;
    end
    imagesc(slice(v));
    colormap(cmap);
    title(title_str);
    if ~isempty(cax)
        caxis(cax);
    end
    xlabel(xlbl);    
    ylabel(ylbl);
    ax = gca;
    ax.XTick = xtick;
    ax.XTickLabel = xticklabel;
    ax.YTick = ytick;
    ax.YTickLabel = yticklabel;
    ax.FontSize = 14;    
    set(gcf,'Units','normalized');
    %set(gca,'LooseInset',get(gca,'TightInset'))
    %set(gcf,'OuterPosition',[0.1 0.1 0.5 0.5]);
    set(gcf,'PaperPositionMode','auto');
    if cbar, colorbar; end
end