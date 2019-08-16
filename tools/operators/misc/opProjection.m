classdef opProjection < opSpot
% OPPROJECTION - Orthogonal projection on to subspaces. Projects a vector or matrix on to span(A)
% or span(A)^{\perp}, given a basis for the subspace, A
%    
% Curt Da Silva
% curtd@math.ubc.ca
%
% Usage:
%    P = opProjection(A, {perp},{isOrthog})
%
% Input:
%    A         - basis for the subspace (or SPOT operator)
%    perp      - true: projects on to span(A)^{\perp}
%                false: projects on to span(A) (default)
%    isOrthog  - true: treats A as orthogonal, if known apriori (skips checking). Leads to much faster projections if true.
%                false: A is not known to be orthogonal and compensates for this fact in the projection (slower)
%
% Output:
%    P         - SPOT operator for the subspace projector
    properties(SetAccess = protected)
        A,perp,isOrthog,AAinv
    end
    methods        
        function op = opProjection(A,perp,isOrthog)                      
            op = op@opSpot('Subspace projection',size(A,1),size(A,1));
            assert(isnumeric(A) || isa(A,'opSpot'), 'A must be a SPOT operator or a matrix');
            assert(size(A,1) >= size(A,2), 'A must be a square or tall matrix');
            if(exist('perp','var')==0)
                perp = false;
            end
            if(exist('isOrthog','var')==0)
                isOrthog = norm(A' * A - eye(size(A,2))) < 1e-8;
            end                        
            if isnumeric(A)
                op.A = opMatrix(A);
            else
               op.A = A; 
            end
            op.A = A;
            op.perp = perp;
            op.sweepflag = true;
            if isOrthog
                op.AAinv = opDirac(size(A,2));
            else
                op.AAinv = opInverse(A' * A);
            end            
        end
        function out = test(op)
            x = randn(op.n,1);
            y = randn(op.n,1);
            
            z = op.multiply(x,1);
            t = y' * z;
            
            z = op.multiply(y,-1);
            s = x' * z;
            
            e1 = abs(s - t)/abs(max(s,t));
                                         
            if~(e1<1e-10); fprintf(2,'opProjection: adjoint test failed, error = %g\n',e1); end
           
            z = op.multiply(x,1);
            e2 = (x - z)' * z;
            if~(abs(e2)<1e-10); fprintf(2,'opProjection: residual not perpendicular to projection, error=%g\n',e2); end
            
            out = (e1<1e-10) & (abs(e2)<1e-10);
               
       end
    end
    methods(Access = protected)
        function y = multiply(op,x,mode)
            y = (op.A)'* x;
            y = op.AAinv * y;
            y = op.A * y;
            if(op.perp)               
                y = x - y;
            end

            z = (op.A)'*y;
            if norm(z) > 1e-4
                z = op.AAinv * z;
                z = op.A * z;
                if op.perp
                    z = y - z;
                end
                y = z;           
            end
        end
    end
    
end
