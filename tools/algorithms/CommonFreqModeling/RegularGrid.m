classdef RegularGrid 
% Regularly spaced Cartesian grid with optional point restriction
%
% Curt Da Silva, 2016
%
    properties (SetAccess = protected)
        x,y,z,idx;
    end
    
    methods
        % G = RegularGrid(xt,yt,zt,idx)
        % 
        % Input:
        %  xt,yt,zt - 1D regularly spaced arrays (set zt = [] for a 2D grid)
        %  idx - 1D indices to indicate the grid points to restrict to (default: 1:numel(xt)*numel(yt)*numel(zt))
        %
        function G = RegularGrid(xt,yt,zt,idx)
            G.x = xt;
            G.y = yt;
            if exist('zt','var')
                G.z = zt;           
            else
                G.z = [];
            end
            if exist('idx','var')==0||isempty(idx),
                idx = [];
            end
            G.idx = idx;                
        end
        % G.numpts() - number of points in the grid
        % 
        function n = numpts(G)            
            n = length(G.idx);
            if n==0
                if isempty(G.z)
                    n = length(G.x)*length(G.y);
                else
                    n = length(G.x)*length(G.y)*length(G.z);
                end
            end
                
        end
        
        % 
        function op = interp_operator(G,xt,yt,zt,index,interp_type,adj_mode)
            if exist('index','var')==0||isempty(index),index = G.idx; else index = G.idx(index); end
            if exist('adj_mode','var')==0,adj_mode = false; end
            if isempty(G.z)
                if adj_mode
                    op = opInterp(interp_type,G.x,xt,G.y,yt);                    
                else
                    op = opInterp(interp_type,xt,G.x,yt,G.y);
                end
            else
                if adj_mode
                    op = opInterp(interp_type,G.x,xt,G.y,yt,G.z,zt);
                else
                    op = opInterp(interp_type,xt,G.x,yt,G.y,zt,G.z);
                end
            end
            if ~isempty(index)&& ( length(index) < G.numpts() )
                R = opRestriction(G.numpts(),index);
                if adj_mode
                    op = op*R';
                else
                    op = R*op;
                end
            end
        end
        
    end
end