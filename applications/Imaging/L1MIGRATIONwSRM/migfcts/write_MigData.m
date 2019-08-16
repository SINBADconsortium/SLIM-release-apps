function write_MigData(options, dm, dm_RTM, dm_deconv, x, info, renewal_info, model_para)
% Syntax:
% write_MigData(options, dm, dm_aj, x, info, renewal_info, model_para)
%
% Note:
% This function is exclusively called by Migration_with_SRM. The only purpose to
% separate it is to make the main driver more human-readable.
%
% Description:
% write output to files according to user's options
%
% Input list:
% options: a struct that contains many other parameters, corresponding to the
%		input of Migration_with_SRM.
% dm: inverted model perturbation
% dm_RTM: cross-correlation reverse time migration result
% dm_deconv: deconvolutional migration result
% x: solution is sparse space. for diagnostic use only
% info: a struct that saves how residule decreases, etc. for diagnostic use only
% renewal_info: a struct that saves all renewal info. for diagnostic use only
% model_para: all parameters used in modelling, for diagnostic use only
%
% Output list: None. Write data to files.
% 
% Author: Ning Tu
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: Feb/14/2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

savesolmat_flag = options.save_solmat;
solmat_name = [options.results_folder, options.solmat_name];
saveevolve_flag = options.save_evolve;
evolve_name = [options.results_folder, options.evolve_name];
preview_name = [options.preview_folder, options.preview_name];

%% save preview file
save(preview_name,'dm','dm_RTM','dm_deconv','options','model_para');
disp('Preview file saved in matlab native format.')

%% save solution file
if savesolmat_flag
    if saveevolve_flag
        load(evolve_name,'saved')
        sol_updates = saved*options.scale_correction;
        clear saved
    else
        sol_updates = [];
    end
    save(solmat_name,'sol_updates','dm','dm_RTM','dm_deconv','x','info','renewal_info','options','model_para')
end

%% save velocity perturbation files
if ~isempty(options.model_pert_out)
	writeDataType([options.results_folder, options.model_pert_out], dm);
end