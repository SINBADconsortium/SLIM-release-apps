% HELMHOLTZ_3D_DERIVATIVE - Computes the derivative of the helmholtz_3d.m
% matrix
% 
% USAGE:
%   [dH ddH] = helmholtz_3d_derivative(v,comp_grid,freqhz,unit,mode)
%
% INPUT:
%   freqhz  - frequency in Hertz
%   model   - velocity model in meters per second or slowness squared; this is
%             a *3D object* with dimension (nx,ny,nz), NOT A VECTOR!
%   h       - vector containing the grid spacing in each direction
%   unit    - String set either to "meters_per_second" or "slowness_squared"
%
% OUTPUT:
%   dH      - Virtually a tensor, but in practice its just a matrix with
%             dimensions (27,nx*ny*nz)
%
% AUTHOR: Rafael Lago
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: August, 2014
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%
%-----------------------------------------------------------------------------
function [dH,ddH] = helmholtz_3d_derivative(v,comp_grid,freqhz,unit,mode)

MMM=1; NMM=2; PMM=3; MNM=4; NNM=5; PNM=6; MPM=7; NPM=8; PPM=9;
MMN=10;NMN=11;PMN=12;MNN=13;NNN=14;PNN=15;MPN=16;NPN=17;PPN=18;
MMP=19;NMP=20;PMP=21;MNP=22;NNP=23;PNP=24;MPP=25;NPP=26;PPP=27;

opts = struct;
if ~exist('mode','var'), opts.mode = opBandStorage.MULT_ONLY;
else opts.mode = mode;
end
%
% Compute PML
%---------------------------------------------------------
nx = comp_grid.nt(1); ny = comp_grid.nt(2); nz = comp_grid.nt(3);
v  = reshape(v,nx,ny,nz);

% Convert velocity to the DERIVATIVE of the wavenumber
% that's what we use anyway;
%--------------------------------------------------------------
if nargout==1
    [~,v] = param2wavenum(v,freqhz,unit);
else
    [~,v,dv] = param2wavenum(v,freqhz,unit);
end

for i=1:nargout    
    Dwn = zeros(1,nx+2,ny+2,nz+2);
    if i==1
        Dwn(1,2:nx+1,2:ny+1,2:nz+1) = reshape(v,1,nx,ny,nz);
    else
        if norm(vec(dv)) < 1e-14
            ddH = opZeros(size(out,2)); 
            break;
        else
            Dwn(1,2:nx+1,2:ny+1,2:nz+1) = reshape(dv,1,nx,ny,nz);
        end
    end
    
    hxyz = comp_grid.dt(1)^2 + comp_grid.dt(2)^2 +comp_grid.dt(3)^2;
    
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
    
    % Initialize dH with zeroes
    out = zeros(27,nx,ny,nz);
    
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
    %  Just mass lumping terms
    %----------------------------------------
    out(MNN,:,:,:) = - wm2*Dwn(1,nxM, nyN, nzN);
    out(PNN,:,:,:) = - wm2*Dwn(1,nxP, nyN, nzN);
    out(NMN,:,:,:) = - wm2*Dwn(1,nxN, nyM, nzN);
    out(NPN,:,:,:) = - wm2*Dwn(1,nxN, nyP, nzN);
    out(NNM,:,:,:) = - wm2*Dwn(1,nxN, nyN, nzM);
    out(NNP,:,:,:) = - wm2*Dwn(1,nxN, nyN, nzP);
    out(NPP,:,:,:) = - wm3*Dwn(1,nxN, nyP, nzP);
    out(NMM,:,:,:) = - wm3*Dwn(1,nxN, nyM, nzM);
    out(NPM,:,:,:) = - wm3*Dwn(1,nxN, nyP, nzM);
    out(NMP,:,:,:) = - wm3*Dwn(1,nxN, nyM, nzP);
    out(PNP,:,:,:) = - wm3*Dwn(1,nxP, nyN, nzP);
    out(MNM,:,:,:) = - wm3*Dwn(1,nxM, nyN, nzM);
    out(PNM,:,:,:) = - wm3*Dwn(1,nxP, nyN, nzM);
    out(MNP,:,:,:) = - wm3*Dwn(1,nxM, nyN, nzP);
    out(PMN,:,:,:) = - wm3*Dwn(1,nxP, nyM, nzN);
    out(MPN,:,:,:) = - wm3*Dwn(1,nxM, nyP, nzN);
    out(PPN,:,:,:) = - wm3*Dwn(1,nxP, nyP, nzN);
    out(MMN,:,:,:) = - wm3*Dwn(1,nxM, nyM, nzN);
    out(PPP,:,:,:) = - wm4*Dwn(1,nxP, nyP, nzP);
    out(MMM,:,:,:) = - wm4*Dwn(1,nxM, nyM, nzM);
    out(PPM,:,:,:) = - wm4*Dwn(1,nxP, nyP, nzM);
    out(MMP,:,:,:) = - wm4*Dwn(1,nxM, nyM, nzP);
    out(PMP,:,:,:) = - wm4*Dwn(1,nxP, nyM, nzP);
    out(MPM,:,:,:) = - wm4*Dwn(1,nxM, nyP, nzM);
    out(PMM,:,:,:) = - wm4*Dwn(1,nxP, nyM, nzM);
    out(MPP,:,:,:) = - wm4*Dwn(1,nxM, nyP, nzP);
    out(NNN,:,:,:) = - (w1 + 3*w2 + (16*w3a*hxyz)/3 + wm1-1)*Dwn(1,nxN,nyN,nzN);
    
    %
    % Set Zero to the Boundaries
    % Is this needed? - Lago
    %----------------------------
    out(MNN, 1, :, :) = 0;
    out(PNN,nx, :, :) = 0;
    out(NMN, :, 1, :) = 0;
    out(NPN, :,ny, :) = 0;
    out(NNM, :, :, 1) = 0;
    out(NNP, :, :,nz) = 0;
    out(MMN, 1, :, :) = 0; out(MMN, :, 1, :) = 0;
    out(PMN,nx, :, :) = 0; out(PMN, :, 1, :) = 0;
    out(MPN, 1, :, :) = 0; out(MPN, :,ny, :) = 0;
    out(PPN,nx, :, :) = 0; out(PPN, :,ny, :) = 0;
    out(MNM, 1, :, :) = 0; out(MNM, :, :, 1) = 0;
    out(PNM,nx, :, :) = 0; out(PNM, :, :, 1) = 0;
    out(MNP, 1, :, :) = 0; out(MNP, :, :,nz) = 0;
    out(PNP,nx, :, :) = 0; out(PNP, :, :,nz) = 0;
    out(NMM, :, 1, :) = 0; out(NMM, :, :, 1) = 0;
    out(NPM, :,ny, :) = 0; out(NPM, :, :, 1) = 0;
    out(NMP, :, 1, :) = 0; out(NMP, :, :,nz) = 0;
    out(NPP, :,ny, :) = 0; out(NPP, :, :,nz) = 0;
    out(MMM, 1, :, :) = 0; out(MMM, :, 1, :) = 0; out(MMM, :, :, 1) = 0;
    out(PMM,nx, :, :) = 0; out(PMM, :, 1, :) = 0; out(PMM, :, :, 1) = 0;
    out(MPM, 1, :, :) = 0; out(MPM, :,ny, :) = 0; out(MPM, :, :, 1) = 0;
    out(PPM,nx, :, :) = 0; out(PPM, :,ny, :) = 0; out(PPM, :, :, 1) = 0;
    out(MPP, 1, :, :) = 0; out(MPP, :,ny, :) = 0; out(MPP, :, :,nz) = 0;
    out(PPP,nx, :, :) = 0; out(PPP, :,ny, :) = 0; out(PPP, :, :,nz) = 0;
    out(MMP, 1, :, :) = 0; out(MMP, :, 1, :) = 0; out(MMP, :, :,nz) = 0;
    out(PMP,nx, :, :) = 0; out(PMP, :, 1, :) = 0; out(PMP, :, :,nz) = 0;
    
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
    out = reshape(out,27,nx*ny*nz);
    
    out = opBandStorage(out,idx,opts);
    if i==1
        dH = out;
    else
        ddH = out; 
    end
end
end
