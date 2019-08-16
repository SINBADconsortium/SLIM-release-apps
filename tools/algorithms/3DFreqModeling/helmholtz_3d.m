% 3D Helmholtz operator, discretized using 27-point stencil with PML
% [Operto et al. 2007. Geophysics 72(5), SM195]
%
% USE:
%   H = helmholtz_3d(model,h,np,freqhz,unit)
%
% INPUT:
%   freqhz  - frequency in Hertz
%   model   - velocity model in meters per second or slowness squared; this is a 
%             *3D object* with dimension (nx,ny,nz), NOT A VECTOR!
%   h       - 3D vector containing the grid spacing in each direction
%   np      - structure containing the number of np points in each direction 
%             np.x, np.y, np.bottom and np.top (the later two stand for the z 
%             direction). For free-surface, set np.top to zero
%   unit    - String set either to "m/s" (meters per second) or "s2/m2" 
%             (slowness squared)
%
% OUTPUT:
%   H       - array with dimensions (27,nx*ny*nz)
%
% AUTHOR: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: January, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
% 
%-------------------------------------------------------------------------------
function [H idx] = helmholtz_3d(v,h,np,freqhz,unit)

   MMM=1; NMM=2; PMM=3; MNM=4; NNM=5; PNM=6; MPM=7; NPM=8; PPM=9;
   MMN=10;NMN=11;PMN=12;MNN=13;NNN=14;PNN=15;MPN=16;NPN=17;PPN=18;
   MMP=19;NMP=20;PMP=21;MNP=22;NNP=23;PNP=24;MPP=25;NPP=26;PPP=27;
   
   %
   % Compute PML
   %---------------------------------------------------------
   [nx,  ny,  nz] = size(v);
   [xix, xiy, xiz] = create_pml(nx, ny, nz, np);
   
   
   % Convert velocity/slowness squared to wavenumber - that's 
   % what we use anyway;
   %--------------------------------------------------------------   
   v = param2wavenum(v,freqhz,unit);
   wn = zeros(1,nx+2,ny+2,nz+2);
   wn(1,2:nx+1,2:ny+1,2:nz+1) = reshape(v,1,nx,ny,nz);

   % We use only the squared values anyway
   %----------------------------------------------
   hx2 = h(1)^2;
   hy2 = h(2)^2;
   hz2 = h(3)^2;
   
   hyz  = hy2 + hz2;
   hxy  = hy2 + hx2;
   hxz  = hx2 + hz2;
   hxyz = hx2 + hy2 + hz2;
   
   %
   % Choose weights for mass lumping
   %----------------------------------------------
   w1  = 1.8395262e-5;
   w2  = (0.890077)/3;
   w3  = (0.1099046)/4;
   wm1 = 0.49649658;
   wm2 = (0.4510125)/6;
   wm3 = (0.052487)/12;
   wm4 = (0.45523e-5)/8;
   w3a = ((w3*3)/(4*hxyz));
   
   % Initialize H with zeroes
   H = zeros(27,nx,ny,nz);
   
   %--------------------------------------------------------------------------------!
   % DANGER - DANGER - DANGER - DANGER - DANGER - DANGER - DANGER - DANGER - DANGER !
   %--------------------------------------------------------------------------------!
   %                                                                                !
   %                              Stencil Generation                                !
   %                          Indexes, indexes EVERYWHERE                           !
   %                                                                                !
   %--------------------------------------------------------------------------------!
   % DANGER - DANGER - DANGER - DANGER - DANGER - DANGER - DANGER - DANGER -  LAGO  !
   %--------------------------------------------------------------------------------!
   
   [~, H_MNN, ~ ] = meshgrid(1:ny, xix(1:nx,1), 1:nz); H_MNN = reshape(H_MNN,1,nx,ny,nz);
   [~, H_PNN, ~ ] = meshgrid(1:ny, xix(1:nx,2), 1:nz); H_PNN = reshape(H_PNN,1,nx,ny,nz);
   [H_NMN, ~, ~ ] = meshgrid(xiy(1:ny,1), 1:nx, 1:nz); H_NMN = reshape(H_NMN,1,nx,ny,nz);
   [H_NPN, ~, ~ ] = meshgrid(xiy(1:ny,2), 1:nx, 1:nz); H_NPN = reshape(H_NPN,1,nx,ny,nz);
   [~, ~, H_NNM ] = meshgrid(1:ny, 1:nx, xiz(1:nz,1)); H_NNM = reshape(H_NNM,1,nx,ny,nz);
   [~, ~, H_NNP ] = meshgrid(1:ny, 1:nx, xiz(1:nz,2)); H_NNP = reshape(H_NNP,1,nx,ny,nz);
   
   % This helps a lot the indeces with the mass lumping
   %-----------------------------------------------------
   nxM = 1:nx;
   nyM = 1:ny;
   nzM = 1:nz;
   nxN = 2:nx+1;
   nyN = 2:ny+1;
   nzN = 2:nz+1;
   nxP = 3:nx+2;
   nyP = 3:ny+2;
   nzP = 3:nz+2;
   
   % 
   %         Main + Stencil X + Stencil Y + Stencil Z         |       Stencil X Mixed        |       Stencil Y Mixed        |      Stencil Z Mixed         |          D stencils           |            Mass          |     D1     |      D2     |     D3      |     D4      |
   %------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
   H(MNN,:,:,:) = - (w1/hx2 + w2/hx2 + w2/hxz + w2/hxy)*H_MNN                                + (w2/(2*hxz))*(H_NNM + H_NNP) + (w2/(2*hxy))*(H_NMN + H_NPN) - 8*w3a*H_MNN                   - wm2*wn(1,nxM, nyN, nzN); %   2*MNN         2*MNN         2*MNN         2*MNN     
   H(PNN,:,:,:) = - (w1/hx2 + w2/hx2 + w2/hxz + w2/hxy)*H_PNN                                + (w2/(2*hxz))*(H_NNM + H_NNP) + (w2/(2*hxy))*(H_NMN + H_NPN) - 8*w3a*H_PNN                   - wm2*wn(1,nxP, nyN, nzN); %   2*PNN         2*PNN         2*PNN         2*PNN     
   H(NMN,:,:,:) = - (w1/hy2 + w2/hy2 + w2/hyz + w2/hxy)*H_NMN + (w2/(2*hyz))*(H_NNM + H_NNP)                                + (w2/(2*hxy))*(H_MNN + H_PNN) - 8*w3a*H_NMN                   - wm2*wn(1,nxN, nyM, nzN); %   2*NMN         2*NMN         2*NMN         2*NMN     
   H(NPN,:,:,:) = - (w1/hy2 + w2/hy2 + w2/hyz + w2/hxy)*H_NPN + (w2/(2*hyz))*(H_NNM + H_NNP)                                + (w2/(2*hxy))*(H_MNN + H_PNN) - 8*w3a*H_NPN                   - wm2*wn(1,nxN, nyP, nzN); %   2*NPN         2*NPN         2*NPN         2*NPN     
   H(NNM,:,:,:) = - (w1/hz2 + w2/hz2 + w2/hxz + w2/hyz)*H_NNM + (w2/(2*hyz))*(H_NMN + H_NPN) + (w2/(2*hxz))*(H_MNN + H_PNN)                                - 8*w3a*H_NNM                   - wm2*wn(1,nxN, nyN, nzM); %   2*NNM         2*NNM         2*NNM         2*NNM     
   H(NNP,:,:,:) = - (w1/hz2 + w2/hz2 + w2/hxz + w2/hyz)*H_NNP + (w2/(2*hyz))*(H_NMN + H_NPN) + (w2/(2*hxz))*(H_MNN + H_PNN)                                - 8*w3a*H_NNP                   - wm2*wn(1,nxN, nyN, nzP); %   2*NNP         2*NNP         2*NNP         2*NNP     
   H(NPP,:,:,:) =                                             - (w2/(2*hyz))*(H_NNP + H_NPN)                                                               + 2*w3a*(H_MNN + H_PNN)         - wm3*wn(1,nxN, nyP, nzP); %     MNN           PNN
   H(NMM,:,:,:) =                                             - (w2/(2*hyz))*(H_NNM + H_NMN)                                                               + 2*w3a*(H_MNN + H_PNN)         - wm3*wn(1,nxN, nyM, nzM); %     MNN           PNN
   H(NPM,:,:,:) =                                             - (w2/(2*hyz))*(H_NPN + H_NNM)                                                               + 2*w3a*(H_MNN + H_PNN)         - wm3*wn(1,nxN, nyP, nzM); %     MNN                                      PNN
   H(NMP,:,:,:) =                                             - (w2/(2*hyz))*(H_NMN + H_NNP)                                                               + 2*w3a*(H_MNN + H_PNN)         - wm3*wn(1,nxN, nyM, nzP); %     MNN                                      PNN
   H(PNP,:,:,:) =                                                                            - (w2/(2*hxz))*(H_NNP + H_PNN)                                + 2*w3a*(H_NMN + H_NPN)         - wm3*wn(1,nxP, nyN, nzP); %     NMN                        NPN
   H(MNM,:,:,:) =                                                                            - (w2/(2*hxz))*(H_NNM + H_MNN)                                + 2*w3a*(H_NMN + H_NPN)         - wm3*wn(1,nxM, nyN, nzM); %     NMN                        NPN
   H(PNM,:,:,:) =                                                                            - (w2/(2*hxz))*(H_NNM + H_PNN)                                + 2*w3a*(H_NMN + H_NPN)         - wm3*wn(1,nxP, nyN, nzM); %                   NMN                        NPN
   H(MNP,:,:,:) =                                                                            - (w2/(2*hxz))*(H_NNP + H_MNN)                                + 2*w3a*(H_NMN + H_NPN)         - wm3*wn(1,nxM, nyN, nzP); %                   NMN                        NPN
   H(PMN,:,:,:) =                                                                                                           - (w2/(2*hxy))*(H_PNN + H_NMN) + 2*w3a*(H_NNM + H_NNP)         - wm3*wn(1,nxP, nyM, nzN); %                                NNM           NNP
   H(MPN,:,:,:) =                                                                                                           - (w2/(2*hxy))*(H_MNN + H_NPN) + 2*w3a*(H_NNM + H_NNP)         - wm3*wn(1,nxM, nyP, nzN); %                                NNM           NNP
   H(PPN,:,:,:) =                                                                                                           - (w2/(2*hxy))*(H_PNN + H_NPN) + 2*w3a*(H_NNM + H_NNP)         - wm3*wn(1,nxP, nyP, nzN); %     NNM           NNP
   H(MMN,:,:,:) =                                                                                                           - (w2/(2*hxy))*(H_MNN + H_NMN) + 2*w3a*(H_NNM + H_NNP)         - wm3*wn(1,nxM, nyM, nzN); %     NNM           NNP
   H(PPP,:,:,:) =                                                                                                                                          - 2*w3a*(H_NPN + H_NNP + H_PNN) - wm4*wn(1,nxP, nyP, nzP); % NPN + NNP      NNP + PNN    NPN + PNN    
   H(MMM,:,:,:) =                                                                                                                                          - 2*w3a*(H_NMN + H_NNM + H_MNN) - wm4*wn(1,nxM, nyM, nzM); % NMN + NNM      NNM + MNN    NMN + MNN    
   H(PPM,:,:,:) =                                                                                                                                          - 2*w3a*(H_NNM + H_PNN + H_NPN) - wm4*wn(1,nxP, nyP, nzM); % NNM + PNN      NNM + PNN                  NPN + PNN
   H(MMP,:,:,:) =                                                                                                                                          - 2*w3a*(H_NNP + H_MNN + H_NMN) - wm4*wn(1,nxM, nyM, nzP); % NNP + MNN      NMN + NNP                  NMN + MNN
   H(PMP,:,:,:) =                                                                                                                                          - 2*w3a*(H_PNN + H_NMN + H_NNP) - wm4*wn(1,nxP, nyM, nzP); % PNN + NMN                   NMN + NNP     PNN + NNP
   H(MPM,:,:,:) =                                                                                                                                          - 2*w3a*(H_MNN + H_NPN + H_NNM) - wm4*wn(1,nxM, nyP, nzM); % MNN + NPN                   NPN + NNM     MNN + NNM
   H(PMM,:,:,:) =                                                                                                                                          - 2*w3a*(H_PNN + H_NMN + H_NNM) - wm4*wn(1,nxP, nyM, nzM); %                PNN + NMN    PNN + NNM     NMN + NNM
   H(MPP,:,:,:) =                                                                                                                                          - 2*w3a*(H_MNN + H_NPN + H_NNP) - wm4*wn(1,nxM, nyP, nzP); %                MNN + NPN    MNN + NNP     NPN + NNP
   
   
   H(NNN,:,:,:) = -(w1 + 3*w2 + (16*w3a*hxyz)/3 + wm1-1)*wn(1,nxN,nyN,nzN) ...
                      + (w1/hx2 + w2/hx2 + w2/hxz + w2/hxy + 8*w3a)*(H_MNN + H_PNN) ...
                      + (w1/hy2 + w2/hy2 + w2/hyz + w2/hxy + 8*w3a)*(H_NMN + H_NPN) ...
                      + (w1/hz2 + w2/hz2 + w2/hyz + w2/hxz + 8*w3a)*(H_NNM + H_NNP);


   % 
   % Set Zero to the Boundaries
   % Is this needed? - Lago
   %-------------------------
   H(MNN, 1, :, :) = 0;
   H(PNN,nx, :, :) = 0;
   H(NMN, :, 1, :) = 0;
   H(NPN, :,ny, :) = 0;
   H(NNM, :, :, 1) = 0;
   H(NNP, :, :,nz) = 0;
   H(MMN, 1, :, :) = 0; H(MMN, :, 1, :) = 0;
   H(PMN,nx, :, :) = 0; H(PMN, :, 1, :) = 0;
   H(MPN, 1, :, :) = 0; H(MPN, :,ny, :) = 0;
   H(PPN,nx, :, :) = 0; H(PPN, :,ny, :) = 0;
   H(MNM, 1, :, :) = 0; H(MNM, :, :, 1) = 0;
   H(PNM,nx, :, :) = 0; H(PNM, :, :, 1) = 0;
   H(MNP, 1, :, :) = 0; H(MNP, :, :,nz) = 0;
   H(PNP,nx, :, :) = 0; H(PNP, :, :,nz) = 0;
   H(NMM, :, 1, :) = 0; H(NMM, :, :, 1) = 0;
   H(NPM, :,ny, :) = 0; H(NPM, :, :, 1) = 0;
   H(NMP, :, 1, :) = 0; H(NMP, :, :,nz) = 0;
   H(NPP, :,ny, :) = 0; H(NPP, :, :,nz) = 0;
   H(MMM, 1, :, :) = 0; H(MMM, :, 1, :) = 0; H(MMM, :, :, 1) = 0;
   H(PMM,nx, :, :) = 0; H(PMM, :, 1, :) = 0; H(PMM, :, :, 1) = 0;
   H(MPM, 1, :, :) = 0; H(MPM, :,ny, :) = 0; H(MPM, :, :, 1) = 0;
   H(PPM,nx, :, :) = 0; H(PPM, :,ny, :) = 0; H(PPM, :, :, 1) = 0;
   H(MPP, 1, :, :) = 0; H(MPP, :,ny, :) = 0; H(MPP, :, :,nz) = 0;
   H(PPP,nx, :, :) = 0; H(PPP, :,ny, :) = 0; H(PPP, :, :,nz) = 0;
   H(MMP, 1, :, :) = 0; H(MMP, :, 1, :) = 0; H(MMP, :, :,nz) = 0;
   H(PMP,nx, :, :) = 0; H(PMP, :, 1, :) = 0; H(PMP, :, :,nz) = 0;
   
   %
   % Write down the IDX vector
   %--------------------------------------------------------------
   idx = zeros(1,27);
   idx([MMP NMP PMP MNP NNP PNP MPP NPP PPP]) = idx([MMP NMP PMP MNP NNP PNP MPP NPP PPP]) + nx*ny;
   idx([MMM NMM PMM MNM NNM PNM MPM NPM PPM]) = idx([MMM NMM PMM MNM NNM PNM MPM NPM PPM]) - nx*ny;
   idx([MPM NPM PPM MPN NPN PPN MPP NPP PPP]) = idx([MPM NPM PPM MPN NPN PPN MPP NPP PPP]) + nx;
   idx([MMM NMM PMM MMN NMN PMN MMP NMP PMP]) = idx([MMM NMM PMM MMN NMN PMN MMP NMP PMP]) - nx;
   idx([PMM PNM PPM PMN PNN PPN PMP PNP PPP]) = idx([PMM PNM PPM PMN PNN PPN PMP PNP PPP]) + 1;
   idx([MMM MNM MPM MMN MNN MPN MMP MNP MPP]) = idx([MMM MNM MPM MMN MNN MPN MMP MNP MPP]) - 1;

   % Vectorize second dimension
   %---------------------------------------------------------------
   H = reshape(H,27,nx*ny*nz);
end 
