classdef DataGrid
    properties (Access = protected)
        src_grid, rec_grid;
    end
    methods
        function D = DataGrid(s_grid,r_grid)
            D.src_grid = s_grid;
            D.rec_grid = r_grid;
        end        
        
        function nsrc = numsrcs(G)
            nsrc = G.src_grid.numpts();
        end
                
        function nrec = numrecs(G,index)
            if isa(G.rec_grid,'cell')
                nrec = G.rec_grid{index}.numpts();
            else
                nrec = G.rec_grid.numpts();
            end
        end
        function [Ps,Pr] = interp_operators(G,index,comp_grid,interp_type)
            if exist('interp_type','var')==0, interp_type = 'sinc'; end
            mode_2d = comp_grid.nt(3)==1; adj_mode = true;
            if mode_2d
                [zt,xt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);
                Ps = G.src_grid.interp_operator(zt,xt,[],[],interp_type,adj_mode);
                if isa(G.rec_grid,'cell')
                    Pr = G.rec_grid{index}.interp_operator(zt,xt,[],[],interp_type,adj_mode)';
                else
                    Pr = G.rec_grid.interp_operator(zt,xt,[],[],interp_type,adj_mode)';
                end
            else
                [xt,yt,zt] = odn2grid(comp_grid.ot,comp_grid.dt,comp_grid.nt);
                Ps = G.src_grid.interp_operator(xt,yt,zt,[],interp_type,adj_mode);
                if isa(G.rec_grid,'cell')
                    Pr = G.rec_grid{index}.interp_operator(xt,yt,zt,[],[],interp_type,adj_mode)';
                else
                    Pr = G.rec_grid.interp_operator(xt,yt,zt,[],interp_type,adj_mode)';
                end
            end

           
        end
    end
end