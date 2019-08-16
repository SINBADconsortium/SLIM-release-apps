function test_suite = test_DoubleKaczmarzSweep
% Author: Art Petrenko
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
%         
% Date: July, 2014

% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

% For background on the Kaczmarz sweeps and CARP, please see the following papers: 
% Bjorck and Elfving, "Accelerated projection methods for computing pseudoinverse solutions of systems of linear equations", 1979
% Gordon and Gordon, "Component-averaged row projections: a robust, block-parallel scheme for sparse linear systems", 2005
% Gordon and Gordon, "CARP-CG: A robust and efficient parallel solver for linear systems, applied to strongly convection dominated PDEs", 2010

	test_suite=buildFunctionHandleTestSuite(localfunctions);
end % function test_DoubleKaczmarzSweep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function freshFixture = setupDefineLinearSystem()
% Setup function that returns a linear system for testing.
%
% freshFixture.R	 - 27 x prod(n) array with non-zero diagonals of the matrix H as rows 
% freshFixture.idx   - 1 x 27 vector of diagonal offsets
% freshFixture.R_7_point_star	  - 7 x prod(n) array with non-zero diagonals of the matrix H as rows 
% freshFixture.idx_7_point_star   - 1 x 27 vector of diagonal offsets
% freshFixture.x     - prod(n) x 1 initial guess
% freshFixture.n	 - dimensions of system
% freshFixture.matrixType - cell array of strings with the names of the matrix types to facilitate access to these matrices by the test subfunctions.

	freshFixture.n = [10,11,12];
	initial_guess_type = 'zeros';

	N = prod(freshFixture.n);

	% Generate matrix for the 27-point cube stencil
	stencil_length = 27;
	freshFixture.R_27_point_cube = randn([3,3,3, freshFixture.n]) + 1i*randn([3,3,3, freshFixture.n]);
	% Set to zero those elements that would be zero in a domain of size n
	freshFixture.R_27_point_cube(1   , :   , :   , 1   , :   , :)   = 0;
	freshFixture.R_27_point_cube(end , :   , :   , end , :   , :)   = 0;
	freshFixture.R_27_point_cube(:   , 1   , :   , :   , 1   , :)   = 0;
	freshFixture.R_27_point_cube(:   , end , :   , :   , end , :)   = 0;
	freshFixture.R_27_point_cube(:   , :   , 1   , :   , :   , 1)   = 0;
	freshFixture.R_27_point_cube(:   , :   , end , :   , :   , end) = 0;
	freshFixture.R_27_point_cube = reshape(freshFixture.R_27_point_cube,[stencil_length,N]);
	row_norm_inverses = 1./sqrt(sum(abs(freshFixture.R_27_point_cube).^2,1));
	freshFixture.R_27_point_cube = bsxfun(@times,freshFixture.R_27_point_cube,row_norm_inverses);
	% Calculate offsets of diagonals of system matrix stored in R
	freshFixture.idx_27_point_cube = zeros(1,stencil_length);
	j = 1;
	for i_slow = -1:1
		for i_med = -1:1
			for i_fast = -1:1
				freshFixture.idx_27_point_cube(j) = i_fast + i_med*freshFixture.n(1) + i_slow*freshFixture.n(1)*freshFixture.n(2);
				j = j + 1;
			end
		end
	end
	freshFixture.matrixType{1} = '27_point_cube';

	% Generate matrix for the 7-point star stencil
	stencil_length = 7;
	freshFixture.R_7_point_star = randn([stencil_length, freshFixture.n]) + 1i*randn([stencil_length, freshFixture.n]);
	% Set to zero those elements that would be zero in a domain of size n.
	freshFixture.R_7_point_star(1   , :   , :   , 1)   = 0;
	freshFixture.R_7_point_star(7   , :   , :   , end) = 0;
	freshFixture.R_7_point_star(2   , :   , 1   , :)   = 0;
	freshFixture.R_7_point_star(6   , :   , end , :)	  = 0;
	freshFixture.R_7_point_star(3   , 1   , :   , :)   = 0;
	freshFixture.R_7_point_star(5   , end , :   , :)   = 0;
	freshFixture.R_7_point_star = reshape(freshFixture.R_7_point_star,[stencil_length,N]);
	row_norm_inverses = 1./sqrt(sum(abs(freshFixture.R_7_point_star).^2,1));
	freshFixture.R_7_point_star = bsxfun(@times,freshFixture.R_7_point_star, row_norm_inverses);
	% Calculate offsets of diagonals of system matrix stored in R
	freshFixture.idx_7_point_star = zeros(1,stencil_length);
	j = 1;
	for i_slow = -1:1
		switch i_slow
		case 0
			for i_med = -1:1
				switch i_med
				case 0
					for i_fast = -1:1
						freshFixture.idx_7_point_star(j) = i_fast + i_med*freshFixture.n(1) + i_slow*freshFixture.n(1)*freshFixture.n(2);
						j = j + 1;
					end
				otherwise
					freshFixture.idx_7_point_star(j) = i_med * freshFixture.n(1);
					j = j + 1;
				end
			end
		otherwise
			freshFixture.idx_7_point_star(j) = i_slow * freshFixture.n(1) * freshFixture.n(2); 
			j = j + 1;
		end
	end
	freshFixture.matrixType{2} = '7_point_star';

	% Initial guess
	if strcmp(initial_guess_type,'randn')
		freshFixture.x = randn(N,1) + 1i*randn(N,1);
	elseif strcmp(initial_guess_type,'zeros')
		temp = 1 + 1i;
		freshFixture.x = zeros(N,1,'like',temp);
	elseif strcmp(initial_guess_type,'ints')
		freshFixture.x = (1:2:2*N).' + 1i*(2:2:2*N).';
	else
		fprintf('Invalid selection for initial guess type!\n');
		freshFixture.x = [];
	end

end % function setupDefineLinearSystem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test case functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_DoubleCARPSweepMatrixHermitian(freshFixture)
% Adjoint test for CARP.
	for i = 1:length(freshFixture.matrixType)
		matrix_type = freshFixture.matrixType{i};
		fprintf('Matrix type %s... ', matrix_type);
		fprintf('Testing that the double Kaczmarz/CARP sweeps are Hermitian...\n');

		N = length(freshFixture.(['R_' matrix_type]));
		b = zeros(N,1,'like',freshFixture.x);

		% The test runs over several parallel modes, as well as the sequential mode
		% (Kaczmarz sweeps). It tests both the 27-point cube stencil and the
		% 7-point star stencil.
		np = 4;
		for ip = 1:np
			fprintf('%d blocks... ', ip);
			% CARP in R^n is equivalent to Kaczmarz in the superspace where the vectors
			% are extended with halos (Gordon and Gordon, 2005). 
			% For sequential Kaczmarz sweeps, the "CGMN system" (Bjorck and
			% Elfving, 1979, equation 4.7) is I - Q, where Q is the double Kaczmarz
			% sweep operator for the original system with a right hand side of
			% zero. 
			IminusQfun = @(x)(last1Halo(x - dsweep(freshFixture.(['R_' matrix_type]), freshFixture.(['idx_' matrix_type]), x, b, 1.0, ip), ip, max(abs(freshFixture.(['idx_' matrix_type])))));

			x = randn(N,1);
			y = randn(N,1);

			assertElementsAlmostEqual(last1Halo(y, ip, max(abs(freshFixture.(['idx_' matrix_type]))))' * IminusQfun(x), ...
									  IminusQfun(y)' * last1Halo(x, ip, max(abs(freshFixture.(['idx_' matrix_type])))), 'relative', 1e-6);
		end
		fprintf('\n');
	end
end

function test_DoubleCARPSweepMatrixPosDef(freshFixture)
% Use the Cholesky decomposition or the full eigenvalue spectrum to test for
% positive definiteness. The Cholesky decomposition is recommneded by
% Mathworks, but the chol function incorrectly flags some matrices as
% non-positive definite if they have small (1e-8) positive eigenvalues. Using
% the eig function to find all the eigenvalues seems to work around this, but
% also takes much longer.
	for i = 1:length(freshFixture.matrixType)
		matrix_type = freshFixture.matrixType{i};
		fprintf('Matrix type %s... ', matrix_type);

		N = length(freshFixture.(['R_' matrix_type]));
		x = zeros(N,1,'like',freshFixture.x);
		b = zeros(N,1,'like',freshFixture.x);

		% Construct matrix representation by operating on columns of the identity
		fprintf('\nTesting that the double Kaczmarz/CARP sweeps are positive semidefinite... \n');
		np = 4;
		for ip = 1:np
			fprintf('%d blocks... ', ip);
			haloGridpointOrder = last1Halo(1:N, ip, max(abs(freshFixture.(['idx_' matrix_type]))))';
			IminusQ = complex(zeros(length(haloGridpointOrder)), zeros(length(haloGridpointOrder)));
			for i = haloGridpointOrder
				x(i) = 1;
				IminusQ(:,i) = last1Halo(x - dsweep(freshFixture.(['R_' matrix_type]), freshFixture.(['idx_' matrix_type]), x, b, 1.0, ip), ip, max(abs(freshFixture.(['idx_' matrix_type]))));
				x(i) = 0;
			end

		%	fprintf('Computing Cholesky decomposition to test positive definiteness.\n');
		%	[~,p] = chol(IminusQ);
		%	assert(p == 0, 'Cholesky decomposition found that the double Kaczmarz sweep function does not correspond to a positive definite matrix, p = %d.', p);

			eigenvalueSpectrum = eig(IminusQ);
			eigenvalueSpectrum = sort(eigenvalueSpectrum, 'ascend');
			assert(eigenvalueSpectrum(1) >= 0, ...
				'The eigenvalue spectrum of the double CARP sweep function (%d blocks) indicates it is not positive definite. Smallest eigenvalue = %g.', ...
				ip, eigenvalueSpectrum(1));
		end
		fprintf('\n');
	end
end

function test_DoubleCARPSweepMatrixErrorReducing(freshFixture)
% Kaczmarz sweep should monotonically reduce the pseudoinverse error ||A^+ b - x_k||_2
	for i = 1:length(freshFixture.matrixType)
		matrix_type = freshFixture.matrixType{i};
		fprintf('Matrix type %s... ', matrix_type);

		fprintf('Testing that the double Kaczmarz sweeps reduce the error... \n');
		N = length(freshFixture.(['R_' matrix_type]));
		n_runs = 10;

		% Original (Helmholtz) system matrix
		A = spdiags((freshFixture.(['R_' matrix_type])).', -(freshFixture.(['idx_' matrix_type])), N, N);
		A = A.';

		np = 4;
		for ip = 1:np
			fprintf('%d blocks... ', ip);
			for i = 1:n_runs
				x_soln = complex(randn(N,1), randn(N,1));
				b = A*x_soln;
				x_DKSWP = dsweep(freshFixture.(['R_' matrix_type]), freshFixture.(['idx_' matrix_type]), freshFixture.x, b, 1.0, ip);
				assert(norm(last1Halo(x_soln - x_DKSWP, ip, max(abs(freshFixture.(['idx_' matrix_type]))))) < norm(last1Halo(x_soln - freshFixture.x, ip, max(abs(freshFixture.(['idx_' matrix_type]))))), ...
					'Double CARP sweep (%d blocks) did not reduce the error.', ip); 
			end
		end
		fprintf('\n');
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x = dsweep(R,idx,x,b,w,np)
% Helper function that runs the double Kaczmarz sweep. Note that the last
% argument to sweepR_mex, np, is optional.
	x = sweepR_mex(R,idx,x,b,w,1,np);
	x = sweepR_mex(R,idx,x,b,w,-1,np);
end

function [xHalos, cumPartitionHalos] = last1Halo(xNoHalos, np, haloWidth)
% Helper function that takes a vector and returns a longer vector with halos
% inserted into it.  
%
% INPUT
% xNoHalos	- original vector.
% np		- the number of blocks (parallel processors) that the domain will be split into.
% haloWidth - the number of grid points to include in the halo. 
% 
% OUTPUT
% xHalos 			- vector with halos inserted.
% cumPartitionHalos - vector with the index of the last element of each block in xHalos. There is an additional first element: 0.
	xNoHalos = xNoHalos(:);
	N = length(xNoHalos);
	xHalos 	 = zeros([N + 2*haloWidth*(np-1), 1], 'like', xNoHalos); 
	partition = [repmat(ceil(N/np), [1, rem(N, np)]), ...
				 repmat(floor(N/np), [1, np - rem(N, np)])];
	cumPartition = [0, cumsum(partition)];
	% The first and last blocks have one halo each, the rest of the blocks have
	% two halos.
	if np == 1
		numHalos = [0];
	else
		numHalos = [1, repmat(2, [1, max(np-2,0)]), 1];
	end
	cumPartitionHalos   = cumPartition + [0, cumsum(numHalos)]*haloWidth;
	for ip = 1:np
		xHalos(cumPartitionHalos(ip)+1 : cumPartitionHalos(ip+1)) = ...
			xNoHalos(cumPartition(ip)+1 - (ip>1)*haloWidth : cumPartition(ip+1) + (ip<np)*haloWidth);
	end
	xHalos = xHalos(:);
end

