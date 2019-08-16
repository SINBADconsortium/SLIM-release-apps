function hss_ind = HSS_cell(A, split_levels)

%--------------------------------------------------------------------------------
% Factorize matrix in HSS structure
% HSS_cell(A, split_levels) will return the split indices of the matrix A for
% the number of split_levels, sorted in order from largest to smallest chunk.
% 
% Use:
%   HSS_cell(A, split_levels)
%
% Input: 
%                A - matrix A   
%     split_levels - number of levels of partitioning 

% Authors: Thomas lai and Rajiv Kumar
%          Seismic Laboratory for Imaging and Modeling
%          Department of Earth, Ocean, and Atmospheric Sciences
%          The University of British Columbia
%  
% Date: 10 May 2013

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%--------------------------------------------------------------------------------

% Variables
[m,n] = size(A);
mat_queue = [1 m 1 n]'; % queue for sub-matrices, start with A
hss_ind = []; % n x 4 matrix for storing indices in the form [ms,me,ns,ne;...  transposed

% For split levels
for i = 1:split_levels

    new_mat_queue = [];

    % Loop through all the sub chunks
    for j = 1:size(mat_queue,2)

        current_chunk = mat_queue(:,j);
        
        % Bottom left chunk index
        msb = ceil((current_chunk(1)+current_chunk(2))/2);
        meb = current_chunk(2);
        nsl = current_chunk(3);
        nel = floor((current_chunk(3)+current_chunk(4)-1)/2);
        hss_ind = [hss_ind [msb;meb;nsl;nel]]; % append to result
        
        % Top right
        mst = current_chunk(1);
        met = floor((current_chunk(1)+current_chunk(2)-1)/2);
        nsr = ceil((current_chunk(3)+current_chunk(4))/2);
        ner = current_chunk(4);
        hss_ind = [hss_ind [mst;met;nsr;ner]]; % append to result
        
        % Process next level of chunks
        % Top left
        new_mat_queue = [new_mat_queue [mst;met;nsl;nel]];
        
        % Bottom right
        new_mat_queue = [new_mat_queue [msb;meb;nsr;ner]];
    end

    mat_queue = new_mat_queue;

end % split_levels for loop

hss_ind = [hss_ind mat_queue]; % append diagonals

end % function

