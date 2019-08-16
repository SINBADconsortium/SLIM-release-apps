function [varargout] = PDEfunc_compute_wri(func,A,Dobs,U,Q,src_weights,Pr,H,lambda,dm,to_phys,idx_m)
% PDEfunc_compute_wri - Routine to perform the PDEfunc computations for WRI 
% 
% Curt Da Silva, 2016
    
    S = [lambda*Q;Dobs];
    G = A'*A;
    if isempty(U)
        U = A\S;
    end
    U = U*src_weights;
    Dp = Pr*U;
    pde_res = H*U-Q;    
    sum_srcs = @(x) to_phys*sum(real(x),2);
    switch func
      case PDEopts.FIELD
        varargout = {U};
      case PDEopts.FORW_MODEL
        varargout = {Dp};
      case PDEopts.OBJ
        fdata = 0.5 * norm((Dp-Dobs),'fro')^2;
        fpde = 0.5* lambda^2*norm(pde_res,'fro')^2;
        f = fdata + fpde;
        if nargout >= 2
            V = dH'*pde_res;
            g = lambda^2 * sum_srcs(T(U)'*V);
            ddH1   = ddH(idx_m,idx_m);
            dH1    = dH(idx_m,idx_m);
            n = size(to_phys,1);
            for i=1:size(U,2)          
                U_diag  = U(:,i); U_diag = spdiags(U_diag(idx_m),0,n,n);
                res_diag = pde_res(:,i); res_diag = spdiags(res_diag(idx_m),0,n,n);
                O   = ddH1'*res_diag;
                P   = dH1*U_diag;                              
                h   = lambda^2 * real( U_diag'*O +  P'*P);
            end
            varargout = {f,g,h};
        else
            varargout = {f};
        end
      case PDEopts.HESS_GN
        dHdm   = dH*opDiag_swp(dm);
        dU     = lambda^2 * ( G \ (-dHdm'*(pde_res) - H' * (dHdm*U)));
        Z1     = Pr*dU; Z2 = lambda * (dHdm * U + H * dU);
        Z1     = Pr' * Z1 + lambda * H' * Z2;
        Z1     = lambda^2 * ( G \ Z1 );
        R      = conj( -dH' * (pde_res)) .* Z1;
        R      = sum(R - conj(U)  .* (dH' * H * Z1),2);
        Z2     = sum( conj( U ) .* (dH' * (lambda*Z2)) , 2);
        varargout =  {to_phys*(real(R + Z2))};
      case PDEopts.HESS
        dHdm   = dH*opDiag_swp(dm);
        ddHtdm  = ddH'*opDiag_swp(dm);
        dU     = lambda^2 * (G \ (-dHdm'*(pde_res) - H' * (dHdm*U)));
        output = lambda^2*to_phys*sum(real( conj(U) .* (dH'*(dHdm*U)) ),2);
        % d/dm d/du phi(m,u) [du]
        output = output + lambda^2*to_phys*sum(real( conj(U) .* (ddHtdm*pde_res)),2);
        output = output + lambda^2*to_phys*sum(real( conj(dU) .* (dH' *pde_res)  ),2);
        output = output + lambda^2*to_phys*sum(real( conj(U) .* ( dH' *(H*dU) ) ),2);   
        varargout = {output};
    end
end