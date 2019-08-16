function [v,model] = load_model(model_path, ns, model_name, entire_model)
% LOAD_MODEL - Load a 3D velocity model
%  Curt Da Silva, 2016
%
% Usage:
%   [v,model] = load_model(model_path,ns,model_name,entire_model);
%
%
if exist('entire_model','var')==0
    entire_model = false;
end

model = struct;
switch model_name
    case 'edam'
        vlow = 2000; vhi = 2200;
        if exist('ns','var')
            nx = ns(1); ny = ns(2); nz = ns(3);
        else
            nx = 100; ny = 100; nz = 100;
        end
        o = [0,0,0];        
        d = [25,25,25];
        [X,Y,Z] = ndgrid(1:nx,1:ny,1:nz);
        R = nx/4;
        v = vlow*ones(nx,ny,nz);
        v( (X-nx/2).^2 + (Y-ny/2).^2 + (Z-ny/2).^2 < R^2 ) = vhi;
        
    case 'compass'
        nx = 1911; ny = 2730; nz = 341;
        o = [0,0,0];
        dx = 25; dy = 25; dz = 12;
        d = [dx,dy,dz];
        
        % read velocity
        v = ReadSegyFast(model_path);
        v = reshape(v,nz,nx,ny);
        v = permute(v,[2 3 1]);                
        
    case {'overthrust','overthrust_small'}
        path = fileparts(model_path); addpath(path);       
        [v,n,d,o] = rsf_read_all(model_path);
        nx = n(1); ny = n(2); nz = n(3);        
        v = reshape(v,nx,ny,nz);
    otherwise
        error('Unrecognized model');
end

if exist('ns','var')==0 || isempty(ns)
    nx_s = nx; ny_s = ny; nz_s = nz;
else
    nx_s = ns(1); ny_s = ns(2); nz_s = ns(3);
end

model.unit  = 'm/s';
ds = d;
if entire_model
    Lx = opLInterp1D(linspace(0,1,nx),linspace(0,1,nx_s));
    Ly = opLInterp1D(linspace(0,1,ny),linspace(0,1,ny_s));
    Lz = opLInterp1D(linspace(0,1,nz),linspace(0,1,nz_s));
    ds = d .* [nx,ny,nz]./ns;
    v = reshape( opKron(Lz,Ly,Lx)*vec(v), [nx_s,ny_s,nz_s]);
else
    if strcmp(model_name,'overthrust')
        Ix = 1:nx_s;
        Iy = 1:ny_s;
        Iz = 1:nz_s;
        v = v(Ix,Iy,Iz);
    elseif strcmp(model_name,'compass') || strcmp(model_name,'overthrust_small')        
        midx = round(nx/2); midy = round(ny/2); 
        Ix = (midx-nx_s/2):(midx+nx_s/2-1);
        Iy = (midy-ny_s/2):(midy+ny_s/2-1);
        Iz = 1:nz_s;
        v = v(Ix,Iy,Iz);                        
    end
end
model.o = o; model.d = ds; model.n = [nx_s,ny_s,nz_s];


end