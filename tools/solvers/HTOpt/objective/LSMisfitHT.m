function [fk,gk] = LSMisfitHT(x, I, b, dimTree)
% LSMISFITHT - Least-squares misfit for interpolating Hierarchical
%              Tucker tensors and its Riemannian gradient
%
% i.e.      1/|\Omega| * \| P_{\Omega} \phi(x) - b\|_2^2
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
    
    %treeidx2linear = @(h,k) 2^(h-1)-1 + k;               
    stop_idx = length(dimTree.nodes{depth});
        
    
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
    %Pad U,Bs with zeros for cache-efficiency
    for i=1:length(B)
        for j=1:length(B{i})
            if ~isempty(B{i}{j})
                if i==1
                    B{i}{j} = padarray(B{i}{j},K*ones(1,2)-size(B{i}{j}));
                    
                else
                    B{i}{j} = padarray(B{i}{j},K*ones(1,3)-size(B{i}{j}));
                end
                if nargout == 2
                    dB{i}{j} = zeros(size(B{i}{j}));
                end
                if i > 1                    
                    Bperm{i}{j} = matricize(B{i}{j},3,[1 2]);
                else
                    Bperm{i}{j} = B{i}{j};
                end
            end
        end
    end
    for i=1:length(U)
        U{i} = padarray(U{i},[0,K-size(U{i},2)]);
        if nargout == 2
            dU{i} = zeros(size(U{i}));
        end
    end    
    %M = sqrt(size(I,2));            
    M = 1;
    if isdistributed(I) || iscodistributed(I)
        %Parallel version
        if nargout == 1
            spmd
                codist = getCodistributor(b);
                Iloc = getLocalPart(I);
                bloc = getLocalPart(b);
                rloc = zeros(length(bloc),1);
                
                Utmp = zeros(nnodes,K); 
                Uidx = zeros(1,d);            
                
                for i=1:length(U)
                    Uidx(i) = 2^(leaf_indices(1,i)-1)-1+leaf_indices(2,i);
                end
                for ii=1:size(Iloc,2)  
                    for i=1:length(U)
                        Utmp(Uidx(i),:) = U{i}(Iloc(i,ii),:);
                    end
                    Utmp = evaluateU(Utmp,B,stop_idx,depth);
                    rloc(ii) = (Utmp(1,1) - bloc(ii))/M;                     
                end
                partSize = codist.Partition;
                codist = codistributor1d(1,partSize,[size(b,1),1]);       
                r = codistributed.build(rloc,codist,'noCommunication');     
            end
        else            
            Uloc = Composite();
            dUloc = Composite();
            dBloc = Composite();
            Bloc = Composite();
            for i=1:parpool_size()
                Uloc{i} = U;
                dUloc{i} = dU;
                Bloc{i} = B;
                dBloc{i} = dB;
            end
            spmd
                codist = getCodistributor(b);
                Iloc = getLocalPart(I);
                bloc = getLocalPart(b);
                rloc = zeros(length(bloc),1);
                
                Utmp = zeros(nnodes,K); 
                dUtmp = zeros(nnodes,K);       
                Uidx = zeros(1,d);                            
                for i=1:length(U)
                    Uidx(i) = 2^(leaf_indices(1,i)-1)-1+leaf_indices(2,i);
                end
                
                for ii=1:size(Iloc,2)                               
                    for i=1:length(Uloc)
                        Utmp(Uidx(i),:) = Uloc{i}(Iloc(i,ii),:);
                    end
                    Utmp = evaluateU(Utmp,Bloc,stop_idx,depth);
                    rloc(ii) = (Utmp(1,1) - bloc(ii))/M;   
                    dUtmp(1,1) = rloc(ii);
                    [dUtmp,dBloc] = evaluatedU(Utmp,dUtmp,Bperm,dBloc,stop_idx,depth);
                    for i=1:length(dUloc)
                        dUloc{i}(Iloc(i,ii),:) = dUloc{i}(Iloc(i,ii),:) + dUtmp(Uidx(i),:);
                    end
                end                
                dx = cell(1,2);
                dx{1} = dUloc; dx{2} = dBloc;
                dx = gop(@htplus,dx,1);
                r = codistributed.build(rloc,codist,'noCommunication');     
            end
            dx = dx{1};        
            dU = dx{1}; dB = dx{2};            
        end           
    else
        %Serial version
        r = zeros(size(I,2),1);          
        Utmp = zeros(nnodes,K);      
        Uidx = zeros(1,d);
        for i=1:length(U)
            Uidx(i) = 2^(leaf_indices(1,i)-1)-1+leaf_indices(2,i);
        end
        
        dUtmp = zeros(nnodes,K);           
        
        if nargout == 1
            for ii=1:size(I,2)                                            
                %Copy the relevant parts of U in to the temporary array                                
                for i=1:length(U)                    
                    Utmp(Uidx(i),:) = U{i}(I(i,ii),:);
                end
                Utmp = evaluateU(Utmp,B,stop_idx,depth);
                r(ii) = (Utmp(1,1) - b(ii))/M;                   
            end        
        else                                   
            for ii=1:size(I,2)                                            
                %Copy the relevant parts of U in to the temporary array                             
                for i=1:length(U)
                    Utmp(Uidx(i),:) = U{i}(I(i,ii),:);
                end
                Utmp = evaluateU(Utmp,B,stop_idx,depth);
                r(ii) = (Utmp(1,1) - b(ii))/M;                 
                
                %Fill up (implicit) derivative tree
                dUtmp(1,1) = r(ii);
                [dUtmp,dB] = evaluatedU(Utmp,dUtmp,Bperm,dB,stop_idx,depth);
                for i=1:length(dU)
                    dU{i}(I(i,ii),:) = dU{i}(I(i,ii),:) + dUtmp(Uidx(i),:);
                end                
            end                       
        end                
    end       
    fk = 0.5 * norm(r,2)^2;
    if nargout == 2 
        %Clip dU,dB back to their original sizes
        for i=1:length(dU)
            r = dimTree.rank(leaf_indices(1,i),leaf_indices(2,i));
            dU{i} = dU{i}(:,1:r);
        end
        numB = 1;
        for i=1:length(B)
            for j=1:length(B{i})
                if ~isempty(B{i}{j}) 
                    k = dimTree.rank(i,j);
                    kl = dimTree.rank(i+1,2*(j-1)+1);
                    kr = dimTree.rank(i+1,2*(j-1)+2);
                    if i==1
                        dB{i}{j} = dB{i}{j}(1:kl,1:kr);
                    else                                        
                        dB{i}{j} = dB{i}{j}(1:kl,1:kr,1:k);
                        numB = numB+1;
                    end
                end
            end
        end            
        gk = dimTree.toVec(dU,dB);
        gk = project_horizontal(x,gk,dimTree);
    end
end

function Utmp = evaluateU(Utmp,B, stop_idx, depth)
    for i=depth-1:-1:1
        if i==depth-1
            endj = stop_idx/2;
        else
            endj = 2^(i-1);
        end
        for j=1:endj
            this_idx = 2^(i-1)-1 + j;
            left_idx = 2^(i)-1 + 2*(j-1)+1;
            right_idx = left_idx+1;

            Btmp = B{i}{j};
            [kl,kr,k] = size(Btmp);
            UL = Utmp(left_idx,:); UR = Utmp(right_idx,:);
            if i == 1
                Uval = (UL * Btmp) * UR';
            else
                Uval = UR*reshape(UL * reshape(Btmp,[kl,kr*k]),[kr,k]);
            end
            Utmp(this_idx,:) = Uval;
        end
    end
end

function [dUtmp, dB] = evaluatedU(Utmp, dUtmp, B, dB, stop_idx, depth)
    numB = 1;
    for i=1:depth-1
        if i==depth-1
            endj = stop_idx/2;
        else
            endj = 2^(i-1);
        end
        for j=1:endj
            this_idx = 2^(i-1)-1 + j;
            left_idx = 2^(i)-1 + 2*(j-1)+1;
            right_idx = left_idx+1;

            Btmp = B{i}{j};
            
            UL = Utmp(left_idx,:); UR = Utmp(right_idx,:);

            if i==1                      
                dUL = (dUtmp(1,1)* UR) * Btmp';
                dUR = (dUtmp(1,1) * UL) * Btmp;
                dB{i}{j} = dB{i}{j} + (dUtmp(1,1)*(UL')) * UR;
            else
                [k,krkl] = size(Btmp);
                kr = sqrt(krkl); kl = sqrt(krkl);
                dUparent = dUtmp(this_idx,:); 
                Y = reshape(dUparent * Btmp,[kl,kr]);
                dUR = UL * Y;
                dUL = UR * Y';                
                dB{i}{j} = dB{i}{j} + reshape(bsxfun(@times,reshape(bsxfun(@times,UL',UR),kl*kr,1) ,dUparent),kl,kr,k); 
                numB = numB+1;
            end
            dUtmp(left_idx,:) = dUL; dUtmp(right_idx,:) = dUR;
        end
    end
end

function Z = htplus(X,Y)
    dU = X{1}; dB = X{2}; dV = Y{1}; dC = Y{2};
    for i=1:length(dU)
        dU{i} = dU{i} + dV{i};
    end    
    for i=1:length(dB)
        for j=1:length(dB{i})
            if ~isempty(dB{i}{j})
                dB{i}{j} = dB{i}{j} + dC{i}{j};
            end
        end
    end
    Z = cell(1,2);    
    Z{1} = dU;
    Z{2} = dB;
end
