function output = Dv_HT(x, dimTree, w, mode)
% compute the Dw or D^H*w with HT factors once you represent data D in HT
% format 
% 
% Usage :
%        output = Dv_HT(x,dimTree,w,mode)
% Input :
%        x       : the vetorized HT factors
%        dimTree : the structure of HT
%        w       : the probing vector
%        mode    :
%                 '1': D^H*w
%                 '2': D*w
% Output :
%        the output of the data matrix in mutiplication with a vector
%
% Author: Yiming Zhang
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%
% Date: March, 2017.


% form small matrices and internal tensors from their verterized format
[U,B] = dimTree.fromVec(x);

% define the size of the dimensions
nrx     = size(U{1},1);
nsx     = size(U{2},1);
nry     = size(U{3},1);
nsy     = size(U{4},1);

switch mode
    % when you want D^H*v
    case '1'
        if size(w,1) ~= nrx*nry
            error('Wrong size of the probing vector, pls check your size');
        else
            x1     = opKron(U{2},U{1})*matricize(B{2}{1},[1 2])*B{1}{1};
            y1     = opKron(U{4},U{3})*matricize(B{2}{2},[1 2]);
            x1     = reshape(x1,size(U{1},1),size(U{2},1),size(B{1}{1},2));

           for i  = 1:nsx
                temp(i,:) = permute(w,[2 1])*conj(reshape(squeeze(x1(:,i,:))*y1.',nrx*nry,nsy));
           end

           output = temp(:);
        end

    % when you want D*v
    case '2'
        if size(w,1) ~= nsx*nsy
             error('Wrong size of the probing vector, pls check your size');
        else
            output = zeros(nrx*nry,1);
            x1     = opKron(U{2},U{1})*matricize(B{2}{1},[1 2])*B{1}{1};
            y1     = opKron(U{4},U{3})*matricize(B{2}{2},[1 2]);
           for j = 1:nsx
                temp     = x1(nrx*(j-1)+1:nrx*j,:)*y1.';
                temp     = reshape(temp,[nrx*nry,nsy])*w(j:nsx:(nsy-1)*nsx+j);
                output   = output + temp;
           end
        end
end

end