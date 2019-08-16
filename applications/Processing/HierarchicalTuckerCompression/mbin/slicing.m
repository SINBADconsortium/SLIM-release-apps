function gather = slicing(dimTree,x,index,geometry,mode)
% Extracts any common shot gathers or receiver gathers from compressed hierarchical Tucker parameters
%
% Usage:
%   gather = slicing(dimTree,x,index,geometry,mode);
%
% Input:
%   x       -    Vectorized HT parameters
%   dimTree -    Generate U ( HT leaf bases) and B (HT interior nodes)
%
%   index   -    The index you want to extract from src or rec location
%                You have two options here:
%
%                I : A vector of the indeces. For example, you have 121
%                shots in total for a acquisition area [0:10m:100m] x
%                [0:10m:100m], and you want to extract the shots at [0m,0m],
%                [10m,20m],[20m,30m],so the input index vector here should
%                be [1,24,36]^T
%
%                II: Directly the locations of the shots/receivers. Take the
%                same example as above, the input should  be [0,0; 10,20; 20,30]
%
%   geometry -   The structure contains the information of acquisition
%   settings for shots/receivers.
%                geometry.so : origin of source coodinate 
%                geometry.ro : origin of receiver coodinate
%                geometry.sd : interval of source grids
%                geometry.rd : interval of receiver grids
%                geometry.sn : no of grids for source 
%                geometry.rn : no of grids for receiver        
%
%   mode    -    1   :Common shot gathers
%                2   :Common receiver gathers
%
% Output:
%   gather  -    the common shot gather or receiver gather in vecterized
%   form
%
% Note: Pls be careful about the coordinates of your original input data.
%       If keep the orignial data D in canonical organization (receiver x,
%       receiver y, source x, source y), the mode you use will match this
%       function, meaning "1" for common shot gathers.
%
% Author: Yiming Zhang
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%
% Date: March, 2017.


if nargin == 1
    error('Not enough inputs'); 
end
if nargin == 2
    error('Not enough inputs'); 
end
if nargin == 3
    error('Not enough inputs'); 
end
if nargin == 4
    error('Not enough inputs'); 
end
if nargin == 5
    
    % form small matrices and internal tensors from their verterized format
    [U,B] = dimTree.fromVec(x); 
end

    
% algorithm of on-the-fly shots/receivers generation from compressed HT parameters  
switch mode
    % if you want to extract shots 
    case '1' 
        
        % compute the true index along x and y direction 
        if size(index,2) == 1
            int      = zeros(size(index,1),2);
            int(:,2) = floor((index -1)/geometry.sn(1)) + 1;
            int(:,1) = mod((index -1),geometry.sn(1)) + 1;
            index    = int;
            clear int 
        elseif size(index,2) == 2
            index(:,1) = (index(:,1) - geometry.so(1))/geometry.sd(1) + 1;
            index(:,2) = (index(:,2) - geometry.so(2))/geometry.sd(2) + 1;
        else
            error('Wrong input for indeces');
        end
        
         x1     = opKron(U{2}(index(:,1),:),U{1})*matricize(B{2}{1},[1 2]);
         y1     = opKron(U{4}(index(:,2),:),U{3})*matricize(B{2}{2},[1 2]);
         % if you only want to extract one common shot gather 
         if size(index,1) == 1
             gather  = x1*B{1}{1}*y1.';
             gather  = gather(:);
         % if you want to extract multiple shots simultanously 
         else
             num   =   size(index,1);
             nrecx =   size(U{1},1);
             nrecy =   size(U{3},1);
             x1    =   reshape(x1,nrecx,num,rank(dimTree,2,1));
             y1    =   reshape(y1,nrecy,num,rank(dimTree,2,2));
             for i = 1:num
                  int       = squeeze(x1(:,i,:))*B{1}{1}*squeeze(y1(:,i,:)).';
                gather(:,i) = int(:);
             end
          clear int 
         end
         
    % if you want extract receivers 
    case '2'
        
        % compute the true index along x and y direction 
        if size(index,2) == 1
            int      = zeros(size(index,1),2);
            int(:,2) = floor((index -1)/geometry.rn(1)) + 1;
            int(:,1) = mod((index -1),geometry.rn(1)) + 1;
            index    = int;
            clear int 
        elseif size(index,2) == 2
            index(:,1) = (index(:,1) - geometry.ro(1))/geometry.rd(1) + 1;
            index(:,2) = (index(:,2) - geometry.ro(2))/geometry.rd(2) + 1;
        else
            error('Wrong input for indeces');
        end
         
        % if you only want to extract one common receiver gather 
          if size(index,1) == 1
             x1     = opKron(U{2},U{1}(index(:,1),:))*matricize(B{2}{1},[1 2]);
             y1     = opKron(U{4},U{3}(index(:,2),:))*matricize(B{2}{2},[1 2]);
             gather = x1*B{1}{1}*y1.';
             gather = gather(:);
         % if you want to extract multiple receivers simultanously 
          else
              num   =   size(index,1);
              nsrcx =   size(U{2},1);
              nsrcy =   size(U{4},1);
              for i = 1:num
                   x1          = opKron(U{2},U{1}(index(i,1),:))*matricize(B{2}{1},[1 2]);
                   y1          = opKron(U{4},U{3}(index(i,2),:))*matricize(B{2}{2},[1 2]);                                         
                   int         = x1*B{1}{1}*y1.';
                   gather(:,i) = int(:);
              end
           clear int    
         end
end
end

