% Plot slices of 3D FWI results
%
% Curt Da Silva, 2016
% 
% Usage:
%  plot_fwi_exp(fwi_experiments,save_dir,model_dir,fig_dir, plots);
%
% Input:
%   fwi_experiments    - vector of fwi_experiment numbers to plot (all corresponding to the same model)
%   save_dir           - directory where results are stored
%   model_dir          - directory where the model is stored
%   fig_dir            - directory where to save figures (if empty, dont save figures)
%   plots              - 2*numplots x 1 cell array, each two entries of which consist of 
%                        'slice_direction', slice_depth, e.g.,
%                        {'x',2500,'z',500,'y',5000}

function plot_fwi_exp(fwi_experiments,save_dir,model_dir,fig_dir,plots)

if exist('fwi_experiments','var')==0
   error('Need fwi_experiments');
end
eval(['fwi_exp' num2str(fwi_experiments(1)) ';']);
v2struct(fwi_params_struct);

switch modelname
  case 'overthrust_small'
    overthrust_small_preamble;    
  case 'edam'
    edam_preamble;
  otherwise
    error('Unrecognized model');
end

for i=1:length(fwi_experiments)
    eval(['fwi_exp' num2str(fwi_experiments(i)) ';']);
    save_file = [save_dir 'fwi_results_' modelname '_exp' num2str(fwi_experiments(i)) '.mat'];    
    results = load(save_file);
    eval(['vexp' num2str(i) ' = results.vhist;']);
end
if isempty(fig_dir), save_fig = @(x) x; 
else save_fig = @(x) print(gcf,'-dpdf',[fig_dir modelname '/' x '.pdf']); end
%% Load model

[v,model] = load_model([model_dir model_file], [nx_s,ny_s,nz_s],modelname,entire_model);
nx = model.n(1); ny = model.n(2); nz = model.n(3);

ns = initialv_smooth;
v0 = reshape(opKron(opSmooth(nz,ns),opSmooth(ny,ns),opSmooth(nx,ns))*vec(v),model.n);


if exist('plots','var')==0 || isempty(plots)        
    switch modelname
      case 'overthrust_small'
        plots = {'z',500,'z',1000,'z',2000,'x',2500,'x',1000,'y',3000};
    end
end
titles = true;

for i=1:length(plots)/2
    figure; 
    nplots = length(fwi_experiments);
    mid = floor(nplots/2);

        
    opts = struct;
    opts.slice_dir = plots{2*i-1};
    opts.slice_depth = plots{2*i};

    subplot(2,nplots,mid);    
    if titles, opts.title_str = 'v true'; end
    vel_plot(v,model,opts); save_fig(['slice_' opts.slice_dir '_' num2str(opts.slice_depth) '_v']);
    ax = caxis;
    opts.cax = ax;
    if titles, opts.title_str = 'v0'; end
    subplot(2,nplots,mid+1);
    vel_plot(v0,model,opts);  save_fig(['slice_' opts.slice_dir '_' num2str(opts.slice_depth) '_v0']);
    for j=1:length(fwi_experiments)
        subplot(2,nplots,nplots+j);
        eval(['last_idx = find(~cellfun(@isempty,vexp' num2str(j) '),1,''last'');']);
        if titles, opts.title_str = ['vest - exp ' num2str(fwi_experiments(j))]; end
        eval(['vel_plot(vexp' num2str(j) '{last_idx},model,opts); save_fig([''slice_'' opts.slice_dir ''_'' num2str(opts.slice_depth) ''_vexp' num2str(fwi_experiments(j)) '_idx_' num2str(last_idx) ''']);']);
    end
    axes('Units','Normal');
    h = title([opts.slice_dir ' slice @ ' num2str(opts.slice_depth) 'm']);
    set(gca,'visible','off');
    set(h,'visible','on');
end

