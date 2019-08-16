function y = HmtimesTest2(op,x,mode)
            G       = op.G; V = op.V; EIG = op.EIG;  AA = op.AA; Px = op.Px;params=op.params; %lambda = op.lambda;
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
                        lambda  = params.lambda(j);
                        % AAloc   = AA{j};
                        EIGi    = EIG(:,j);
                        Vloc    = V(:,:,j);
                        EIGi    = (conj(EIGi) .* EIGi) ./ (conj(EIGi) .* EIGi + lambda^2);
                        SIG     = spdiags(EIGi,0,length(EIGi),length(EIGi));
%                        y       = y + sum(conj(G(:,:,j)) .* (AAloc' *(Vloc' * (SIG * (Vloc * (AAloc * (xx * G(:,:,j))))))), 2);
                        y       = y + lambda^2 * sum(conj(G(:,:,j)) .* (Vloc' * (SIG * (Vloc * (xx * G(:,:,j))))), 2);
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
                for i = 1:length(params);
                    paramst{i} = params{i};
                end
                params = paramst;
                tic;
                k   = 1;
                GAll = [];
                for i = 1:length(V)
                    VV   = V{i};
                    GG   = G{i};
                    EEIG = EIG{i};
                    paramsi = params{i};
                    lambdas = paramsi.lambda;
%                    AAA  = AA{i};

                    
                    for j = 1:size(GG,3)
                        lambda = lambdas(j);
%                        AAloc   = AAA{j};
                        EIGi    = EEIG(:,j);
                        Vloc    = VV(:,:,j);
                        EIGi    = (conj(EIGi) .* EIGi) ./ (conj(EIGi) .* EIGi + lambda^2);
                        SIG     = spdiags(EIGi,0,length(EIGi),length(EIGi));
                        GGG     = GG(:,:,j);keyboard
                        y       = y + lambda^2 * sum(conj(GGG) .* (Vloc' * (SIG * (Vloc * (xx * GGG)))), 2);
                        
                        I_ij{k}    = opNull(size(GG,3));
                        SIG_ij{k}  = opKron(I_ij, SIG);
                        V_ij{k}    = opKron(I_ij, Vloc);
                        for isrc = 1 : size(GGG,2)
                            GAll       = [GAll; opDiag(GGG(:,isrc))];
                        end
                        
                        k       = k + 1;

                    end 
                    
                end
                toc
                
                
                y = real(Px' * y);
                
                
            end
