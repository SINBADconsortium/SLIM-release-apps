function options = Set_Paras(varargin)
% function options = Set_Paras(varargin)
% 
% This function set up struct data.
% Usage:
% 	Suppose a matlab function needs a struct contains parameters a, b, c, d
%   Just do
% 	opts  = Set_Paras('a',1,'b',2,'c',3,'d',4);
% 
% 	if you want to make change to parameters a, b without
% 	change c, d.
% 
% 	opts1 = Set_Paras('a',1,'b',2);
% 	opts  = Set_Paras(opts,opts1);
%
%
% Author: Xiang Li
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: 02, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

i = 1;


while i <= nargin
    arg = varargin{i};
    if ischar(arg)
        
        % A finite state machine to parse name-value pairs.
        if rem(nargin,2) ~= 0
            error('Each input options must be a string/value pair.');
        end
        % get the names of variables
        Names = varargin(1:2:end);
        names = lower(Names);
        options = [];
        m     = length(Names);
        
        for j = 1:m
            eval(['options.' cell2mat(Names(j)) ' = [];']);
        end
        break;
    end
    if ~isa(arg,'struct')
        error(sprintf(['Expected argument %d to be a string parameter name ' ...
            'or an options structure\ncreated with OPTIMSET.'], i));
    end
    if i == 1
        options = arg;
        Names   = struct2nv(options);
    else
        Names1  = struct2nv(arg);
        for j = 1:length(Names1)
            if any(strcmp(Names,deblank(Names1(j,:))))
                eval(['options.',cell2mat(Names1(j,:)),'=  arg.',cell2mat(Names1(j,:)),';']);
            else
                % disp(sprintf('New parameter added ''%s''.', cell2mat(Names1(j,:))));
                Names(end+1) = Names1(j,:);
                eval(['options.',cell2mat(Names1(j,:)),'=  arg.',cell2mat(Names1(j,:)),';']);
            end
        end
    end
    i = i + 1;
end

expval = 0;
while i <= nargin
    arg = varargin{i};
    if ~expval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string parameter name.', i));
        end
        
        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' deblank(Names(j(1),:))];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names(k,:))];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expval = 1;
        
    else
        eval(['options.' cell2mat(Names(j)) '= arg;']);
        expval = 0;
        
    end
    i = i + 1;
end

if expval
    error(sprintf('Expected value for parameter ''%s''.', arg));
end

