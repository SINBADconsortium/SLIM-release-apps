function clear_run_func( x, input, func, save_dir)
%CLEAR_RUN_FUNC - Deletes the stored results from save_run_func.m.
%  
%  Author: Curt Da Silva
% 
%  Usage:
%    clear_run_func( x, input, func, save_dir);
% 
%  Input:
%    x        - vector/matrix with which the function output is associated
%               (usually a helmholtz matrix)
%    
%    input    - single (numeric) input to the function or cell array of input
%               values
% 
%    func     - function handle to an expensive function. If empty, deletes
%               the entire history associated to the array x.
%  
%    save_dir - directory to save results in. Can be empty, but that kind
%               of defeats the whole point of this function.
%     
%

if ~isempty(input) && ~iscell(input), input = {input}; end

if isempty(x), x = input; end

filename = hash_value(x);

filename = filename(1:20);

var_name = hash_value(input);
var_name = genvarname([func2str(func) var_name]);
        
filename = [filename(1:2:end) '_' var_name(1:2:end)];

if isempty(save_dir)
    % Don't do anything
    return;
else       
    % Try to load results of the function        
    if strcmp(save_dir(end),'/')==0, save_dir = [save_dir '/']; end   
    save_file = [save_dir filename '.mat'];
    
    if isempty(input) || isempty(func)
        if exist(save_file,'file')==2
           delete(save_file); 
        else
            files = dir([save_dir filename(1:2:end) '_*.mat']);            
            for i=1:length(files)
                delete([save_dir files{i}.name]);
            end
        end
    end        
   
    % Delete struct in file if found
    if exist(save_file,'file')==2, 
        vars = whos('-file',save_file);
        vname = {vars.name};
        idx_name = find(cellfun(@(x) strcmp(x,var_name),vname)==1);        
        if length(vname)==1,
            delete(save_file);
        elseif ~isempty(idx_name), 
            %Overwrite existing entry with empty value
            var_name = [];
            save(save_file,var_name,'-v7.3','-append');
        end
    end  
end

end

