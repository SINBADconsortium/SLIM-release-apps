function str = matrix2table(A_mean,A_std, col_labels, row_labels, sigfigs)
%MATRIX2TABLE - Converts an input matrix to a table with specified labels.
%  
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
% 
% Usage:
%    str = matrix2table(A_mean,A_std,col_labels,row_labels,sigfigs);
%
% Input:
%    A_mean     - primary values to display in the table (m x n matrix)
%    A_std      - secondary values to display in the table (m x n matrix, can
%                 be empty)
%    col_labels - length n+1 cell array of labels for the columns
%    row_labels - length m cell array of labels for the rows
%    sigfigs    - number of significant digits to display (default: 3)
% 
% Output:
%    str        - cell array of strings, each entry corresponding to one
%                 line of the table

if exist('sigfigs','var')==0
   sigfigs = 3; 
end

if ~isempty(A_std)    
    if size(A_mean,1) ~= size(A_std,1) || size(A_mean,1) ~= size(A_std,1) 
        error('A_mean and A_std must be the same size');
    end        
end
if size(A_mean,1) ~= length(row_labels) || size(A_mean,2)+1 ~= length(col_labels)
    error('Row/column labels must be of the correct size');
end

formatspec = ['%3.' num2str(sigfigs) 'g'];

%Determine maximum column widths
mean_str = arrayfun(@(x) num2str(x,formatspec),A_mean,'UniformOutput',false);
col_width = max(cellfun(@(x) length(x),mean_str),[],1);

if ~isempty(A_std)
    std_str = arrayfun(@(x) num2str(x,formatspec),A_std,'UniformOutput',false);
    col_width = col_width + max(cellfun(@(x) length(x),std_str),[],1) + 3; %compensate for brackets + space
end

heading_widths = cellfun(@(x) length(x),col_labels);
heading_widths(1) = max(heading_widths(1),max(cellfun(@(x) length(x), row_labels)));
col_width = [heading_widths(1),max(col_width, heading_widths(2:end))];

for i=1:size(A_mean,1)
    for j=1:size(A_mean,2)
        if isempty(A_std)
            mean_str{i,j} = [' ' pad_str(mean_str{i,j},col_width(j+1),true) ' '];
        else
            mean_str{i,j} = [' ' pad_str([mean_str{i,j} ' (' std_str{i,j} ')'],col_width(j+1),true) ' '];
        end
    end
end

for i=1:length(col_labels)
    if i < length(col_labels)
        col_labels{i} = [' ' pad_str(col_labels{i},col_width(i),true) ' |'];        
    else
        col_labels{i} = [' ' pad_str(col_labels{i},col_width(i),true) ' '];
    end
end
for i=1:length(row_labels)
    row_labels{i} = [' ' pad_str(row_labels{i},col_width(1),true) ' |'];
end



%Print table

for i=1:size(A_mean,1)+1
    if i==1
        str{i} = [];
        for j=1:length(col_labels)
            if j < size(col_labels(i))
                str{i} = [str{i} col_labels{j} ' |'];
            else
                str{i} = [str{i} col_labels{j}];
            end
        end
    else
        str{i} = row_labels{i-1} ;
        for j=1:size(A_mean,2)
            str{i} = [str{i} mean_str{i-1,j}];
            if j < size(A_mean,2)
                str{i} = [str{i} '|'];
            end
        end
    end
    
end

end

function z = pad_str(y,len,centered)
    if exist('centered','var')==0
        centered = false;
    end
    if centered
        left = ceil((len - length(y))/2);
        right = len - length(y) - left;
        z = [repmat(sprintf(' '),1,left), y, repmat(sprintf(' '),1,right)];
    else
        z = [y, repmat(sprintf(' '),1,len - length(y))];
    end
end