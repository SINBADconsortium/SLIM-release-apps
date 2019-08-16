function [Ps,Pr] = src_rec_interp(ndims,model,comp_grid,src_interp,rec_interp)
% src_rec_interp - Source/Receiver interpolation operators
%
% Curt Da Silva, 2016
% 
% Usage:
%   [Ps,Pr] = src_rec_interp(ndims,model,comp_grid,src_interp,rec_interp);
% 
% Input:
%   ndims      - number of dimensions (2 or 3)
%   comp_grid  - struct returned from discretize_helmholtz
%   src_interp - one of PDEopts.SRC_INTERP_LIN or PDEopts.SRC_INTERP_SINC
%   rec_interp - one of PDEopts.REC_INTERP_LIN or PDEopts.REC_INTERP_SINC
% 
% Output:
%   Ps         - source grid -> computational grid interpolation operator
%   Pr         - computational grid -> receiver grid interpolation operator
%
    SRC_INTERP_LIN = PDEopts.SRC_INTERP_LIN; SRC_INTERP_SINC = PDEopts.SRC_INTERP_SINC;
    REC_INTERP_LIN = PDEopts.REC_INTERP_LIN; REC_INTERP_SINC = PDEopts.REC_INTERP_SINC;

    if ndims==2
        % 2D model is organized as (z,x) for legacy/visualization reasons
        [zt,xt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);           
        if strcmp(src_interp,SRC_INTERP_LIN)
            Ps = opInterp('cubic',zt,model.zsrc,xt,model.xsrc)';            
        elseif strcmp(src_interp,SRC_INTERP_SINC)
            Ps = opInterp('sinc',model.zsrc,zt,model.xsrc,xt);
        end
        if strcmp(rec_interp,REC_INTERP_LIN)
            Pr = opInterp('cubic',zt,model.zrec,xt,model.xrec);
        elseif strcmp(rec_interp,REC_INTERP_SINC)
            Pr = opInterp('sinc',model.zrec,zt,model.xrec,xt)'; 
        end
    else
        % 3D model is (x,y,z)
        [xt,yt,zt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);
        if strcmp(src_interp,SRC_INTERP_LIN)
            Ps = opInterp('cubic',xt,model.xsrc,yt,model.ysrc,zt,model.zsrc)';     
        elseif strcmp(src_interp,SRC_INTERP_SINC)
            Ps = opInterp('sinc',model.xsrc,xt,model.ysrc,yt,model.zsrc,zt);     
        end
        if strcmp(rec_interp,REC_INTERP_LIN)
            Pr = opInterp('cubic',zt,model.zrec,xt,model.xrec);            
        elseif strcmp(rec_interp,REC_INTERP_SINC)
            Pr = opInterp('sinc',model.xrec,xt,model.yrec,yt,model.zrec,zt)';
        end        
    end           
end