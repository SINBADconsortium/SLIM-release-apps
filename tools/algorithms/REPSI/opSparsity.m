function [OP_SPARSITY, emptySolutionVector] = opSparsity(nt, nr, ns, options)
% Define the sparsity frame, which is a 2D Curvelet frame in reciever & shot plane, with wavelets on time axis
        
        CURV = opCurvelet(nr, ns, max(1,ceil(log2(ns) - 3)), 16, 1, options.curvType, 0);
	    WAVE = opSplineWavelet(nt, 1, nt, 1.5, 5);
	    
	    % utility operator that cleans up small imaginary parts of synthesis due to computationally 
	    % inexact thresholding of Curvelet coefficients
        REMOVEIMAG = opEnsureReal(nr*ns);
        CURV = CURV * REMOVEIMAG;
        
    	if options.parallel
    	    OP_SPARSITY = oppKron2Lo(WAVE,CURV); % distributed data is sorted as time-slice gathers
    	    
    	    % initialize a distributed solution vector
    	    emptySolutionVector = dzeros(OP_SPARSITY');
    	    
    	else
    	    OP_SPARSITY = opKron(CURV,WAVE);  % serial data is sorted as shot gathers
    	    
    	    % initialize a local solution vector (empty also works for SPGL1)
    	    emptySolutionVector = [];
    	end
