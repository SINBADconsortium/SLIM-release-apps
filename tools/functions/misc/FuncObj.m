classdef FuncObj < handle
% FuncObj - Anonymous function evaluation without the memory reference penalty
% of anonymous functions. Enables partial specialization of functions at no extra cost/workspace memory references.
%
% Curt Da Silva, 2016
%
% Usage:
%   obj = FuncObj(function_handle,args);
%   [...] = obj(someargs);
% 
% Input:
%   function_handle - handle to function of interest
%   args - n x 1 cell array of arguments, set to [] if the arguments will be specified upon calling obj()
%   someargs - cell array of arguments that fill in the empty arguments in args 
%  
% Output:
%   
        
    properties (SetAccess = protected)
        func,args,arg_labels,idx_empty;
    end
    
    methods 
        function op = FuncObj(function_handle,args,arg_labels)
            op.func = function_handle;
            op.args = args;
            if exist('arg_labels','var'), op.arg_labels = arg_labels; 
            else op.arg_labels = []; end
            idx = find(cellfun(@isempty,args));
            op.idx_empty = idx;
        end
        function varargout = subsref(op,s)
            if length(s)==1
                switch s.type
                  case '.'
                    switch s.subs
                      case 'func', varargout = {op.func};
                      case 'args', varargout = {op.args};
                      case 'arg_labels',varargout = {op.arg_labels};
                      case 'idx_empty',varargout = {op.idx_empty};
                    end
                  case '()'
                    % Actually evaluate the function with the input parameters
                    arguments = op.args;
                    arguments(op.idx_empty) = s.subs(:);
                    varargout{:} = op.func(arguments{:});
                end
            elseif length(s)==2
                switch s(1).type
                  case '.'
                    switch s(1).subs
                      case 'setarg'
                        setarg(op,s(2).subs{:});
                      case 'rmarg'
                        rmarg(op,s(2).subs{1}); 
                      case 'hasarg'
                        varargout{:} = hasarg(op,s(2).subs{1});
                    end
                end
            end
        end
        function setarg(op,i,a)
            if ~isempty(op.arg_labels) && isa(i,'char')
                j = find(strcmp(op.arg_labels,i));
                if ~isempty(j),i=j; end
            end
            op.args{i} = a;
            op.idx_empty = setdiff(op.idx_empty,i);
        end
        function rmarg(op,i)
            op.args{i} = [];
            op.idx_empty = union(op.idx_empty,i);
        end
        function y = hasarg(op,i)
            if isa(i,'char')
                y = ~isempty(cellfun(@(x)strcmp(x,i),op.arg_labels));
            elseif isa(i,'double')
                y = i <= length(op.arg_labels);
            end
        end
    end
end