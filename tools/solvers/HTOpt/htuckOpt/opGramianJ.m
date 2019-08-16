classdef opGramianJ < opSpot
% OPGRAMIANJ - Derivative of the Gramian function (i.e. dimTree.gramian(x))
% Includes the root node, so produces a 0 in that slot for the derivative
%
% Curt Da Silva
% HTOpt v0.1
% curtd@math.ubc.ca
%
% Usage:
%   gramJ = opGramianJ(dimTree,x,{G})
%
% Input:
%   dimTree - dimension tree object
%   x       - current (orthogonalized) point
%   G       - cell array of Gramian matrices corresponding to x, if computed elsewhere (optional)   
        
    properties(SetAccess = protected)
        dimTree,x,G;
    end
    
    methods
        function op = opGramianJ(dimTree,x,G)            
            if(exist('G','var')==0)
                G = dimTree.gramian(x);
            end                  
            op = op@opSpot('Hierarchical Tucker Gramian Jacobian', ...
                           length(dimTree.toVecGram(G)) , length(x));
            op.dimTree = dimTree;
            op.x = x;
           
            op.G = G;            
        end        
    end
    
    methods( Access = protected )
        function y = multiply(op,x,mode)
            dimTree = op.dimTree;
            t = op.x;
            G = op.G;
            [U,B] = dimTree.fromVec(t);
            T = dimTree.emptyTree(U,B);
            T = dimTree.copyCell(T, G, 'G');            
            
            if mode ==1
                [dU,dB] = dimTree.fromVec(x);
                T = dimTree.copyParams(T,dU,dB,'dU','dB');                
                %Forward mode, dx -> dG                                                
                itr = dimTree.iterator('interior_down',T);
               
                sym = @(A) A + A';
                while itr.advance()
                    if itr.isRoot()
                        dG = 0;
                        itr.setValue('dG',0);
                    else
                        dG = itr.getValue('dG');
                    end
                    B = itr.getValue('B'); G = itr.getValue('G');
                    dB = itr.getValue('dB'); 
                    
                    T = ttm(B,dG,3);
                    
                    dGl = ttt(T,B,[2 3],[2 3]) + sym(ttt(ttm(B,G,3),dB,[2 3],[2 3]));
                    dGr = ttt(T,B,[1 3],[1 3]) + sym(ttt(ttm(B,G,3),dB,[1 3],[1 3]));
                    
                    itr.setLeftValue('dG',dGl);
                    itr.setRightValue('dG',dGr);                             
                end
                dG = dimTree.extractParam(itr.T,'dG');
                y = dimTree.toVecGram(dG);
                                               
            else                
                dG = dimTree.fromVecGram(x);                
                T = dimTree.copyCell(T,dG,'dG');
                T = dimTree.copyCell(T,dG,'dGlower');                
                
                itr = dimTree.iterator('interior_up',T);
                while itr.advance()
                   B = itr.getValue('B'); G = itr.getValue('G');
                   dG = itr.getValue('dG');                   
                   dGl = itr.getLeftValue('dGlower'); dGr = itr.getRightValue('dGlower');
                   BG = ttm(B,G,3);
                   dB = ttm(BG,dGl + dGl',1) + ttm(BG,dGr + dGr',2);                   
                                      
                   itr.setValue('dB',dB);                   
                   if ~itr.isRoot()
                       dGlowerL = itr.getLeftValue('dGlower');
                       dGlowerR = itr.getRightValue('dGlower');
                       dGlower = dG + ttt(ttm(B,dGlowerL,1),B,[1 2],[1 2]) + ttt(ttm(B,dGlowerR,2),B,[1 2],[1 2]);                       
                       itr.setValue('dGlower',dGlower);
                   end                                                         
                end
                [dU,dB] = dimTree.extractParams(itr.T,'U','dB');
                for i=1:length(dU)
                   dU{i} = zeros(size(dU{i})); 
                end
                y = dimTree.toVec(dU,dB);
            end
        end
    end
    
    
end