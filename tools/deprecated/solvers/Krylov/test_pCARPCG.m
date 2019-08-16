function test_suite = test_pCARPCG
% Author: Art Petrenko
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: August, 2014

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% For background on the Kaczmarz sweeps and CARP, please see the following papers: 
% Bjorck and Elfving, "Accelerated projection methods for computing pseudoinverse solutions of systems of linear equations", 1979
% Gordon and Gordon, "Component-averaged row projections: a robust, block-parallel scheme for sparse linear systems", 2005
% Gordon and Gordon, "CARP-CG: A robust and efficient parallel solver for linear systems, applied to strongly convection dominated PDEs", 2010

	initTestSuite;
end % function test_pCARPCG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function freshFixture = setupDefineLinearSystem
% Setup function that returns a linear system for testing.
%
% freshFixture.R	 - 27 x prod(n) array with non-zero diagonals of the matrix H as rows (distributed)
% freshFixture.idx   - 1 x 27 vector of diagonal offsets
% freshFixture.x     - prod(n) x 1 initial guess (distributed)
% freshFixture.b     - prod(n) x 1 right-hand side (distributed)
% freshFixture.n 	 - grid dimensions
%
% freshFixture.RGathered, .xGathered and .bGathered are the non-distributed
% copies of the variables.

	n = [10,11,40];
	freshFixture.n = n;
	matrix_type = 'randn';
	initial_guess_type = 'zeros';
	right_hand_side_type = 'randn';
	freshFixture.options.tol = 1e-6;
	freshFixture.options.maxit = 5000;
	freshFixture.options.w     = 1.5;

	N = prod(freshFixture.n);
	stencil_size = [3,3,3];
	stencil_length = prod(stencil_size);

	% System matrix
	if strcmp(matrix_type,'randn')
		spmd
			partitionOfSlowDimension = codistributor1d.defaultPartition(n(end));
			R = randn([stencil_size, n(1:end-1), partitionOfSlowDimension(labindex)]) ... 
		   + 1i*randn([stencil_size, n(1:end-1), partitionOfSlowDimension(labindex)]);
			% Set to zero those elements that would be zero in a domain of size n
			R(1   , :   , :   , 1   , :   , :)   = 0;
			R(end , :   , :   , end , :   , :)   = 0;
			R(:   , 1   , :   , :   , 1   , :)   = 0;
			R(:   , end , :   , :   , end , :)   = 0;
			if labindex == 1
				R(:   , :   , 1   , :   , :   , 1)   = 0;
			end
			if labindex == numlabs
				R(:   , :   , end , :   , :   , end) = 0;
			end
			R = reshape(R,[stencil_length, prod(n(1:end-1))*partitionOfSlowDimension(labindex)]);
			codistributorOfMatrix = codistributor1d(2, prod(n(1:end-1))*partitionOfSlowDimension, [stencil_length,N]);
			Rcodistributed = codistributed.build(R, codistributorOfMatrix);
		end
		freshFixture.R = Rcodistributed;
		freshFixture.RGathered = gather(freshFixture.R);
		% Calculate offsets of diagonals of system matrix stored in R
		freshFixture.idx = zeros(1,stencil_length);
		j = 1;
		for i_slow = -1:1
			for i_med = -1:1
				for i_fast = -1:1
					freshFixture.idx(j) = i_fast + i_med*freshFixture.n(1) + i_slow*freshFixture.n(1)*freshFixture.n(2);
					j = j + 1;
				end
			end
		end
	else
		fprintf('Invalid selection of matrix type!\n');
		freshFixture.R = [];
		freshFixture.idx = [];
	end

	% Initial guess
	spmd
		codistributorOfVector = codistributor1d(1, prod(n(1:end-1))*partitionOfSlowDimension, [N,1]);
	end
	if strcmp(initial_guess_type,'randn')
		spmd
			xCodistributed = randn(N, 1, codistributorOfVector) + 1i*randn(N, 1, codistributorOfVector);
		end
	elseif strcmp(initial_guess_type,'zeros')
		spmd
			xCodistributed = zeros(N, 1, codistributorOfVector) + 1i*zeros(N, 1, codistributorOfVector);
		end
	elseif strcmp(initial_guess_type,'ints')
		spmd
			x = (1 : 2 : 2*codistributorOfVector.Partition(labindex) ).' + 1i*(2 : 2 : 2*codistributorOfVector.Partition(labindex) ).';
			xCodistributed = codistributed.build(x, codistributorOfVector);
		end
	else
		fprintf('Invalid selection for initial guess type!\n');
		freshFixture.x = [];
	end
	freshFixture.x = xCodistributed;
	freshFixture.xGathered = gather(freshFixture.x);

	% Right-hand side
	spmd
		codistributorOfVector = codistributor1d(1, prod(n(1:end-1))*partitionOfSlowDimension, [N,1]);
	end
	if strcmp(initial_guess_type,'randn')
		spmd
			bCodistributed = randn(N, 1, codistributorOfVector) + 1i*randn(N, 1, codistributorOfVector);
		end
	elseif strcmp(initial_guess_type,'zeros')
		spmd
			bCodistributed = zeros(N, 1, codistributorOfVector) + 1i*zeros(N, 1, codistributorOfVector);
		end
	elseif strcmp(initial_guess_type,'ints')
		spmd
			b = (1 : 2 : 2*codistributorOfVector.Partition(labindex) ).' + 1i*(2 : 2 : 2*codistributorOfVector.Partition(labindex) ).';
			bCodistributed = codistributed.build(b, codistributorOfVector);
		end
	else
		fprintf('Invalid selection for right-hand side type!\n');
		freshFixture.b = [];
	end
	freshFixture.b = bCodistributed;
	freshFixture.bGathered = gather(freshFixture.b);

end % function setupDefineLinearSystem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test case functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_solveLinearSystem(freshFixture)
% Verify with a direct solver that pCARPCG can solve linear system distributed
% over several MATLAB workers.

	fprintf('\nTesting that pCARPCG.m (CARP-CG algorithm with parallelization over MATLAB workers) solves a linear system... \n');
	N = prod(freshFixture.n);
	spmd
		nMatlabWorkers = numlabs;
	end
	nMatlabWorkers = nMatlabWorkers{1};

	% Original (Helmholtz) system matrix
	A = spdiags((freshFixture.RGathered).', -(freshFixture.idx), N, N);
	A = A.';
	x_soln = A \ freshFixture.bGathered;

	np = 4;
	fprintf('%d main blocks, ', nMatlabWorkers); 
	for ip = 1:np
		fprintf('%d sub-blocks... ', ip);
		freshFixture.options.n_threads = ip;
		x_pCARPCG = pCARPCG(freshFixture.R,freshFixture.idx,freshFixture.b,freshFixture.x,freshFixture.n,freshFixture.options);
		x_pCARPCG = gather(x_pCARPCG);
		assertVectorsAlmostEqual(x_pCARPCG, x_soln, 'relative', 1e-6);
	end
	fprintf('\n');

end

