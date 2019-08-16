classdef oppJpenApp < oppSpot
% pSPOT wrapper for the penalty method Jacobian. 
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
        mt,D,Q,model,params,G,V,EIG,lambda,AA,Px,distribute,nsrc,nrec,nfreq;
    end
    
    methods
        function op = oppJpenApp(mt,D,Q,model,params,G,V,EIG,lambda,AA,Px,distribute)
           msize1 = 0;
           nfreq  = 0;
           nsrc   = 0;
           nrec   = 0;
           for i = 1:length(G)
               nfreq  = nfreq + size(G{i},3);
            %   nsrc   = nsrc + size(G{i},2);
            %   nrec   = nrec + size(EIG{i},1);
           end
           nrec   = size(EIG{1},1);
           nsrc   = nsrc + size(G{1},2);
           msize1 = nrec * nsrc * nfreq;
           modeltmp = model{1};
           op = op@oppSpot('oppJpenApp',msize1,prod(modeltmp.n));

           op.cflag     = 0;
           op.linear    = 1;
           op.children  = [];
           op.sweepflag = 0;
           op.mt        = mt;
           op.Q         = Q;
           op.D         = D;
           op.model     = model;
           op.params    = params;
           op.G         = G;
           op.V         = V;
           op.EIG       = EIG;
           op.lambda    = lambda;
           op.AA        = AA;
           op.Px        = Px;
           op.distribute = distribute;
           op.nfreq      = nfreq;
           op.nsrc       = nsrc;
           op.nrec       = nrec;
           
           
        end
        
    end
    
    methods ( Access = protected )
        function y = multiply(op,x,mode)
            G       = op.G; V = op.V; EIG = op.EIG;  AA = op.AA; Px = op.Px;params=op.params; opn = op.n; model = op.model;%lambda = op.lambda;
            nfreq   = op.nfreq; nsrc = op.nsrc; nrec = op.nrec;
            
            if mode == 1
                x       = Px * x;
                xx      = spdiags(x, 0, length(x), length(x));
                if op.distribute > 0
                    spmd
                        
                        y      = [];
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
                            EIGi    = sqrt((conj(EIGi) .* EIGi) ./ (conj(EIGi) .* EIGi + lambda^2));
                            SIG     = spdiags(EIGi,0,length(EIGi),length(EIGi));
                            y       = [y; lambda * vec(((SIG * (Vloc * (xx * G(:,:,j)))))) / sigma_j];
                        end 
                    
                        AAloc = []; EIGi = []; Vloc = []; EIGi = []; SIG = []; Gloc = [];
                    end
                    
                    yy  = [];
                    for i = 1:length(y)
                        yy = [yy; y{i}];
                    end
                    y = yy;
                else
                    y      = []; 
                    k      = 1;
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
                            GGG     = GG(:,:,j);
                            y       = [y; lambda * vec(((SIG * (Vloc * (xx * GGG))))) / sigma_j];
                            k       = k + 1;
                        end 
                    end
                
                
                end
            else
                x = reshape(x,nrec,nsrc,nfreq);
                if op.distribute > 0
                    spmd
                        xloc      = x(:,:,params.i_freq);
                        y         = zeros(size(Px,1),1);
                        i         = 1;
                    
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
                            EIGi    = sqrt((conj(EIGi) .* EIGi) ./ (conj(EIGi) .* EIGi + lambda^2));
                            SIG     = spdiags(EIGi,0,length(EIGi),length(EIGi));
                            y       = y + lambda * sum(conj(G(:,:,j)) .* (Vloc' * (SIG * xloc(:,:,j))), 2) / sigma_j;
                        end 
                    
                        y = pSPOT.utils.global_sum(y);
                        AAloc = []; EIGi = []; Vloc = []; EIGi = []; SIG = []; Gloc = [];
                    end
                    y = y{1};
                    y = real(Px' * y);
                else
                    y      = zeros(size(Px,1),1); 
                    k      = 1;
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
                            GGG     = GG(:,:,j);
                            y       = y + lambda * sum(conj(GGG) .* (Vloc' * (SIG * x(:,:,k))), 2) / sigma_j;
                            k       = k + 1;

                        end 
                    end
                    y = real(Px' * y);
                
                
                end
                
                
                
            end
            
        end
    end
end

