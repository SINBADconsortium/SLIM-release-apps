function [fk,gk] = LSMisfitHTmex(x, I, b, dimTree)  
% LSMISFITHTMEX - Least-squares misfit for interpolating Hierarchical
%                 Tucker tensors and its Riemannian gradient, using
%                 a mex implementation.
%
% i.e.       \| P_{\Omega} \phi(x) - b\|_2^2
%
% Input:
%   x       - HT parameter vector
%   I       - d x |\Omega| vector of indices of known data
%             points. For parallel mode, I must be distributed over
%             the second dimension.
%   b       - |\Omega| x 1 vector of data points corresponding to Omega
%   dimTree - dimension tree object    
%
% Output:
%   fk      - objective value
%   gk      - Riemannian gradient
    [U,B] = dimTree.fromVec(x);                        
    [U,Bperm] = dimTree.fromVec(x);
    if nargout == 2
        [dU,dB] = dimTree.fromVec(zeros(length(x),1));
    end
    d = length(dimTree.dims);
    depth = dimTree.depth;    
    
    leaf_indices = dimTree.leaf_indices;                
    
    stop_idx = length(dimTree.nodes{depth-1});
    
    K = 1;
    nnodes = 2^(dimTree.depth)+2*mod(d,2^(dimTree.depth-1))-1;    
    for i=1:length(U)
        K = max(K,size(U{i},2));        
    end
    numB = 0;
    for i=1:length(B)
        for j=1:length(B{i})
            if ~isempty(B{i}{j})
                K = max(K,size(B{i}{j},3));                
                numB = numB + 1;                
            end
        end
    end                  
    Bflat = cell(numB,1);
    idxB = 1;
    for i=1:length(B)
        for j=1:length(B{i})
            if ~isempty(B{i}{j})                                           
                Bflat{idxB} = B{i}{j};
                idxB = idxB + 1;
            end
        end
    end
    for i=1:length(U)      
        U{i} = U{i}.';
        if nargout == 2
            dU{i} = zeros(size(U{i}));
        end
    end             
    
    if isdistributed(I) || iscodistributed(I)
        %Parallel version
        if nargout == 1
            spmd                
                Iloc = getLocalPart(I);
                bloc = getLocalPart(b);
                
                Uidx = zeros(1,d);            
                %Extra -1 for 0-based indexing        
                for i=1:length(U)
                    Uidx(i) = 2^(leaf_indices(1,i)-1)-1+leaf_indices(2,i) - 1;
                end
                Uidx = uint32(Uidx);
                %Serial version                   
                
                fk = mxlsmisfitht(U,Bflat,Iloc,bloc,Uidx,stop_idx);
                %fk = fk * size(Iloc,2)/size(I,2);
                fk = pSPOT.utils.global_sum(fk);    
            end
            fk = fk{1};
        else            
            spmd
                Iloc = getLocalPart(I);
                bloc = getLocalPart(b);                
    
                Uidx = zeros(1,d);            
                %Extra -1 for 0-based indexing        
                for i=1:length(U)
                    Uidx(i) = 2^(leaf_indices(1,i)-1)-1+leaf_indices(2,i) - 1;
                end
                Uidx = uint32(Uidx);
                
                [fk,dUloc,dBloc] = mxlsmisfitht(U,Bflat,Iloc,bloc,Uidx,stop_idx);
                %scaling = size(Iloc,2)/size(I,2);
                scaling = 1;
                fk = fk * scaling;
                for i=1:length(dUloc)
                    dUloc{i} = dUloc{i} * sqrt(scaling);
                end
                for i=1:length(dBloc)
                    dBloc{i} = dBloc{i} * sqrt(scaling);
                end
                fk = pSPOT.utils.global_sum(fk);    
                            
                dx = cell(1,2);
                dx{1} = dUloc; dx{2} = dBloc;                
                dx = gop(@htplus,dx,1);                 
            end
            dx = dx{1};        
            dU = dx{1}; dBin = dx{2};   
            fk = fk{1};
        end           
    else %Serial version            
        Uidx = zeros(1,d);            
        %Extra -1 for 0-based indexing        
        for i=1:length(U)
            Uidx(i) = 2^(leaf_indices(1,i)-1)-1+leaf_indices(2,i) - 1;
        end
        Uidx = uint32(Uidx);        
                       
        if nargout == 1
            fk = mxlsmisfitht(U,Bflat,I,b,Uidx,stop_idx);
        else                                   
            [fk,dU,dBin] = mxlsmisfitht(U,Bflat,I,b,Uidx,stop_idx);
        end                
    end       
    if nargout == 2 
        %Clip dU,dB back to their original sizes
        for i=1:length(dU)
            r = dimTree.rank(leaf_indices(1,i),leaf_indices(2,i));
            dU{i} = dU{i}(1:r,:).';
        end
        
        idxB = 1;
        for i=1:length(B)
            for j=1:length(B{i})
                if ~isempty(B{i}{j}) 
                    k = dimTree.rank(i,j);
                    kl = dimTree.rank(i+1,2*(j-1)+1);
                    kr = dimTree.rank(i+1,2*(j-1)+2);
                    if i==1
                        dB{i}{j} = dBin{idxB}(1:kl,1:kr);
                    else                                        
                        dB{i}{j} = dBin{idxB}(1:kl,1:kr,1:k);                        
                    end
                    idxB = idxB+1;
                end
            end
        end            
        gk = dimTree.toVec(dU,dB);
        gk = project_horizontal(x,gk,dimTree);
    end
end



function Z = htplus(X,Y)
    dU = X{1}; dB = X{2}; dV = Y{1}; dC = Y{2};
    for i=1:length(dU)
        dU{i} = dU{i} + dV{i};
    end    
    for i=1:length(dB)        
        dB{i} = dB{i} + dC{i};
    end
    Z = cell(1,2);    
    Z{1} = dU;
    Z{2} = dB;
end
