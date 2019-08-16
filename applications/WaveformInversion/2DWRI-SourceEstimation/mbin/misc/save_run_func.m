function [ varargout ] = save_run_func( x, input, func, save_dir)
%SAVE_RUN_FUNC - Runs an expensive, deterministic function on an input + saves the
%  result for a future use. If this function is called with identical
%  inputs + function + associated vector, it will return the loaded results
%  from a previous run. Works from within an spmd block.
%  
%  The results from this function can be removed using clear_run_func.m.
%
%  Author: Curt Da Silva
% 
%  Usage:
%    [out1,out2,...] = save_run_func( x, input, func, save_dir);
% 
%  Input:
%    x        - vector/matrix with which the function output is associated
%               (usually a helmholtz matrix)
%    
%    input    - single (numeric) input to the function or cell array of input
%               values
% 
%    func     - function handle to an expensive function
%  
%    save_dir - directory to save results in. Can be empty, but that kind
%               of defeats the whole point of this function.
%    
%  Output:
%    out1,out2,... - outputs of the function
%
%

nout = nargout; 
if ~iscell(input), input = {input}; end

if isempty(x), x = input; end
   
if isempty(save_dir)
    % Don't save output of the function
    output = cell(1,nout);
    [output{:}] = func(input{:});
else       
    filename = hash_value(x);
    filename = filename(1:20);

    var_name = hash_value(input); 
    var_name = genvarname([func2str(func) var_name]);
    
    filename = [filename(1:2:end) '_' var_name(1:2:end)];
    % Try to load results of the function        
    if strcmp(save_dir(end),'/')==0, save_dir = [save_dir '/']; end   
       
    save_file = [save_dir filename '.mat'];          
    
    % Check to see if the old results are in the file and if they're not
    % empty -> don't have to recompute function.
    compute_result = true; 
    if exist(save_file,'file')==2, 
       compute_result = false; 
    end
    
    if compute_result
        output = cell(1,nout);
        [output{:}] = func(input{:});
        % Pack output results in to a struct + save the struct
        eval([var_name ' = struct;']);
        for i=1:nout,             
            eval([ var_name '.( strcat(''a'',num2str(i)) ) = output{i};']);
        end
                       
        if exist(save_file,'file')==2
            save(save_file,var_name,'-v7.3','-append');      
        else
            save(save_file,var_name,'-v7.3');      
        end               
    else 
        y = load(save_file,var_name); y = y.(var_name);        
        output = cell(1,nout);
        for i=1:nout,
           output{i} = y.(['a' num2str(i) ]);
        end 
    end   
end
varargout = output;

end



