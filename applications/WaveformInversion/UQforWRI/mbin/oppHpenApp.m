classdef oppHpenApp < oppSpot
% pSPOT wrapper for the penalty method Hessian. 
% See H_pen.m for more details.
%
% Usage:
%   H = oppHpenApp(m,Q,D,model,params);
%
% Author: Zhilong Fang
% 
% March 2015
%
%
% Author: Zhilong Fang
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: March, 2015
% 
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%
    
    properties 
        mt,D,Q,model,params,G,V,EIG,lambda,AA,Px,distribute;
    end
    
    methods
        function op = oppHpenApp(m,Q,D,model,params)
           op = op@oppSpot('oppHpenApp',length(m),length(m));
           if isfield(params{1}, 'distribute')
               params1 = params{1};
               op.distribute = params1.distribute;
           else
               op.distribute = 1;
           end
           op.cflag     = 0;
           op.linear    = 1;
           op.children  = [];
           op.sweepflag = 0;
           op.mt        = m;
           op.Q         = Q;
           if isa(Q,'Composite') || isa(model,'Composite')
               
           else
               nsrc         = size(Q,2); nrec = length(model.xrec); nfreq = length(model.freq);
               if size(D,2) ~= nfreq
                   D = pSPOT.utils.distVec2distArray(pSPOT.utils.distVectorize(D),[nrec*nsrc,nfreq]); 
               end
           end
           op.D      = D;
           op.model  = model;
           op.params = params;
           if isfield(params{1}, 'Hvel')
                params_1 = params{1};
                Hvel    = params_1.Hvel;
           else
                Hvel = 0;
           end
           
           spmd,
        %       [fStart,fEnd]  = globalIndices(D,2);
               model_loc      = model;
        %       model_loc.freq = model_loc.freq(fStart:fEnd);
        %       if size(Q,3) > 1
        %           Q = getLocalPart(Q);
        %       end  
        %       Dloc             = getLocalPart(D);                                
                Dloc = D; 
                if Hvel > 0
                        m                              = 1./sqrt(m) * 1000;
                        [V G AA EIG Px]  = GenHessMatrixVel(m,Q,Dloc,model_loc,params);
                else
                        [V G AA EIG Px]  = GenHessMatrix(m,Q,Dloc,model_loc,params);
                end
           end
           
           spmd
               for j = 1:size(G,3)
                   V(:,:,j) = V(:,:,j) * AA{j};
               end
           end
           
           op.AA = [];
           
           if op.distribute > 0
               op.V      = V;
               op.G      = G;
               op.EIG    = EIG;
%               op.AA     = AA;
               op.Px     = Px{1};
               % op.lambda = params.lambda;
               op.params = params;
               op.model  = model;
           else
               for i = 1:length(V)
                   VV{i} = V{i};
               end
               op.V   = VV;
               clear VV V
               
               for i = 1:length(G)
                   GG{i} = G{i};
               end
               op.G   = GG;
               clear GG G
               
               for i = 1:length(EIG)
                   EEIG{i} = EIG{i};
               end
               op.EIG   = EEIG;
               clear EIG EEIG
               
               % for i = 1:length(AA)
               %     AAA{i} = AA{i};
               % end
               % op.AA   = AAA;
               % clear AAA AA
               
               op.Px     = Px{1};
               % op.lambda = params.lambda;
               for i = 1:length(params)
                   paramst{i} = params{i};
               end
               op.params = paramst;
               
               for i = 1:length(model)
                   modelt{i} = model{i};
               end
               op.model = modelt;
           end
           
        end
        
        function h = Diags(op) % Obtain the diagonal of the approximated hessian
            G = op.G; V = op.V; EIG = op.EIG; AA = op.AA; Px = op.Px;params = op.params; model = op.model; % ; lambda = op.lambda
            if op.distribute > 0
                spmd
                    h      = zeros(size(op,1),1);
                    i      = 1;
                    for j = 1:size(G,3)
%                    AAloc   = AA{j};
                        if isfield(model,'sigma')
                            sigma_j = model.sigma(j);
                        else
                            sigma_j = 1;
                        end
                        lambda  = params.lambda(j);
                        EIGi    = EIG(:,j);
                        Vloc    = V(:,:,j);
                        EIGi    = sqrt((conj(EIGi) .* EIGi) ./ (conj(EIGi) .* EIGi + lambda^2));
                        SIG     = spdiags(EIGi,0,length(EIGi),length(EIGi));
                        C       = (Vloc' * SIG);
                        
                        for k = 1 : size(G,2)
                            Gloc   = spdiags(G(:,k,j),0,size(G,1),size(G,1));
                            B      = Px' * (Gloc' * C);
                            h      = h + lambda^2 *sum(conj(B).*B,2) / sigma_j^2;
                        
                            % for i = 1: length(h)
                            %     h(i) = h(i) + norm(B(i,:))^2;
                            % end
                        end
                    end 
                    h = pSPOT.utils.global_sum(h);
                    AAloc = []; EIGi = []; Vloc = []; EIGi = []; SIG = []; Gloc = [];
                end
                h = h{1};
                
            else
                h      = zeros(size(op,1),1);
            
                for i = 1:length(V)
                    VV   = V{i};
                    GG   = G{i};
                    EEIG = EIG{i};
                    paramsi = params{i};
                    lambdas = paramsi.lambda;
                    modeli  = model{i};
                    for j = 1:size(GG,3)
                        if isfield(modeli,'sigma')
                            sigma_j = modeli.sigma(j);
                        else
                            sigma_j = 1;
                        end
                        lambda  = lambdas(j);
                        EIGi    = EEIG(:,j);
                        Vloc    = VV(:,:,j);
                        EIGi    = sqrt((conj(EIGi) .* EIGi) ./ (conj(EIGi) .* EIGi + lambda^2));
                        SIG     = spdiags(EIGi,0,length(EIGi),length(EIGi));
                        C       = (Vloc' * SIG);

                        for k = 1:size(GG,2)
                            Gloc   = spdiags(GG(:,k,j),0,size(GG,1),size(GG,1));
                            B      = Px' * (Gloc' * C);
                            h      = h + lambda^2 * sum(conj(B).*B,2) / sigma_j^2;

                        end 
                    end
                
                end
            end
        end
        
        function h = Double(op) % Obtain the diagonal of the approximated hessian
            G = op.G; V = op.V; EIG = op.EIG; AA = op.AA; Px = op.Px;params = op.params; model = op.model; % ; lambda = op.lambda
            if op.distribute > 0
                spmd
                    h      = zeros(size(op,1));
                    i      = 1;
                    for j = 1:size(G,3)
                        if isfield(model,'sigma')
                            sigma_j = model.sigma(j);
                        else
                            sigma_j = 1;
                        end
%                    AAloc   = AA{j};
                        lambda  = params.lambda(j);
                        EIGi    = EIG(:,j);
                        Vloc    = V(:,:,j);
                        EIGi    = sqrt((conj(EIGi) .* EIGi) ./ (conj(EIGi) .* EIGi + lambda^2));
                        SIG     = spdiags(EIGi,0,length(EIGi),length(EIGi));
                        C       = (Vloc' * SIG);
                        
                        for k = 1 : size(G,2)
                            Gloc   = spdiags(G(:,k,j),0,size(G,1),size(G,1));
                            B      = Px' * (Gloc' * C);
                            h      = h + lambda^2 *B*B' / sigma_j^2;
                        
                            % for i = 1: length(h)
                            %     h(i) = h(i) + norm(B(i,:))^2;
                            % end
                        end
                    end 
                    h = pSPOT.utils.global_sum(h);
                    AAloc = []; EIGi = []; Vloc = []; EIGi = []; SIG = []; Gloc = [];
                end
                h = h{1};
                
            else
                h      = zeros(size(op,1));
            
                for i = 1:length(V)
                    VV   = V{i};
                    GG   = G{i};
                    EEIG = EIG{i};
                    paramsi = params{i};
                    lambdas = paramsi.lambda;
                    modeli  = model{i};
                    for j = 1:size(GG,3)
                        if isfield(modeli,'sigma')
                            sigma_j = modeli.sigma(j);
                        else
                            sigma_j = 1;
                        end
                        lambda  = lambdas(j);
                        EIGi    = EEIG(:,j);
                        Vloc    = VV(:,:,j);
                        EIGi    = sqrt((conj(EIGi) .* EIGi) ./ (conj(EIGi) .* EIGi + lambda^2));
                        SIG     = spdiags(EIGi,0,length(EIGi),length(EIGi));
                        C       = (Vloc' * SIG);

                        for k = 1:size(GG,2)
                            Gloc   = spdiags(GG(:,k,j),0,size(GG,1),size(GG,1));
                            B      = Px' * (Gloc' * C);
                            h      = h + lambda^2 * B*B' / sigma_j^2;

                        end 
                    end
                
                end
            end
        end
    end
    
    methods ( Access = protected )
        function y = multiply(op,x,mode)
            G       = op.G; V = op.V; EIG = op.EIG;  AA = op.AA; Px = op.Px;params=op.params; model = op.model;%lambda = op.lambda;
            x       = Px * x;
            xx      = spdiags(x, 0, length(x), length(x));
            if op.distribute > 0
                spmd
                    % if isfield(params,'Threads')
                    %     [maxNumCompThreads(params.Threads)]
                    % end
                    % V = V+1-1;
                    % G = G+1-1;
                    % EIG = EIG + 1 -1;
                    y      = zeros(size(x));
                    i      = 1;
                    
                    for j = 1:size(G,3)
                        if isfield(model,'sigma')
                            sigma_j = model.sigma(j);
                        else
                            sigma_j = 1;
                        end
                        lambda  = params.lambda(j);
                        % AAloc   = AA{j};
                        EIGi    = EIG(:,j);
                        Vloc    = V(:,:,j);
                        EIGi    = (conj(EIGi) .* EIGi) ./ (conj(EIGi) .* EIGi + lambda^2);
                        SIG     = spdiags(EIGi,0,length(EIGi),length(EIGi));
%                        y       = y + sum(conj(G(:,:,j)) .* (AAloc' *(Vloc' * (SIG * (Vloc * (AAloc * (xx * G(:,:,j))))))), 2);
                        y       = y + lambda^2 * sum(conj(G(:,:,j)) .* (Vloc' * (SIG * (Vloc * (xx * G(:,:,j))))), 2) / sigma_j^2;
%                        y       = y + sum(conj(G(:,:,j)) .* (Vloc' * (SIG * ((xx * G(:,:,j))' * Vloc')')), 2);
                        % for k = 1 : size(G,2)
%                             Gloc   = AAloc * spdiags(G(:,k,j),0,size(G,1),size(G,1));
%                             y      = y + lambda^2 * Gloc' * (Vloc' * (SIG * (Vloc * (Gloc * x))));
%                         end
                    end 
                    
                    y = pSPOT.utils.global_sum(y);
                    AAloc = []; EIGi = []; Vloc = []; EIGi = []; SIG = []; Gloc = [];
                end
                
                y = y{1};
                y = real(Px' * y);
%                y = (Px' * y);
            else
                y      = zeros(size(x)); 
                               
                for i = 1:length(V)
                    VV   = V{i};
                    GG   = G{i};
                    EEIG = EIG{i};
                    paramsi = params{i};
                    lambdas = paramsi.lambda;
                    modeli  = model{i};
%                    AAA  = AA{i};
                    for j = 1:size(GG,3)
                        if isfield(modeli,'sigma')
                            sigma_j = modeli.sigma(j);
                        else
                            sigma_j = 1;
                        end
                        lambda = lambdas(j);
%                        AAloc   = AAA{j};
                        EIGi    = EEIG(:,j);
                        Vloc    = VV(:,:,j);
                        EIGi    = (conj(EIGi) .* EIGi) ./ (conj(EIGi) .* EIGi + lambda^2);
                        SIG     = spdiags(EIGi,0,length(EIGi),length(EIGi));
                        GGG     = GG(:,:,j);
%                        y       = y + sum(conj(GG(:,:,j)) .* (Vloc' * (SIG * (Vloc * (xx * GG(:,:,j))))), 2);toc
                        y       = y + lambda^2 * sum(conj(GGG) .* (Vloc' * (SIG * (Vloc * (xx * GGG)))), 2) / sigma_j^2;
%                        y       = y + sum(conj(GGG) .* (Vloc' * (SIG * ((xx * GGG)' * Vloc')')), 2);
%                        y       = y + sum(conj(GG(:,:,j)) .* (AAloc' *(Vloc' * (SIG * (Vloc * (AAloc * (xx * GG(:,:,j))))))), 2);
                        
                        
%                        for k = 1 : size(GG,2)
%
%
%                            Gloc   = AAloc * spdiags(GG(:,k,j),0,size(GG,1),size(GG,1));
%                            y      = y + lambda^2 * Gloc' * (Vloc' * (SIG * (Vloc * (Gloc * x))));
% %                            Gloc   = GG(:,k,j);
% %                            y      = y + lambda^2 * conj(Gloc) .* (AAloc * (Vloc' * (SIG * (Vloc * AAloc * ((Gloc .* x))))));
%
%
%
%                        end
                    end 
                end
                
                y = real(Px' * y);
                
                
            end
            
        end
    end
end

