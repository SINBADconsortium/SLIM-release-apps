function D = F_wri( m, Q, Dobs, model,params )
%F_WRI Forward modelled data for WRI. Useful for visualizing the results of
%   the WRI approach to waveform inversion. Not recommended to be used for
%   computing misfits/whatnot in an inversion algorithm, use misfit_pen or
%   misfit_setup for that. Similar interface to F.m for FWI.
%
% Curt Da Silva, May 2015
% curtd@math.ubc.ca
%
% Usage:
%   D = F_wri( m, Q, Dobs, model, params );
%
% Input:
%   m       - vector with gridded squared slowness in [km^2/s^2]
%   Q       - source matrix. size(Q,1) must match source grid
%             definition, size(Q,2) determines the number of
%             sources, if size(Q,3)>1, it represents a
%             frequency-dependent source and has to be distributed over the last dimension.
%   Dobs    - observed data
%   model   - model struct, see misfit_pen
%   params  - parameter struct, see misfit_pen
%
% Output:
%   D       - Data cube (nrec x nsrc x nfreq) as (distributed) vector. nsrc  = size(Q,2);
%                                                                      nrec  = length(zrec)*length(xrec) 
%                                                                      nfreq = length(freq)
%
    params.wri = true;
    if parpool_size()==0
        D = PDEfunc(PDEopts.FORW_MODEL,m,Q,[],Dobs,model,params);    
    else
        D = PDEfunc_dist(PDEopts.FORW_MODEL,m,Q,[],Dobs,model,params);
    end
D = vec(D);
end

