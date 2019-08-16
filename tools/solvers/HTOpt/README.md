-------------------------------------------------------------------------------
## Hierarchical Tucker Optimization Toolbox v1.0
-------------------------------------------------------------------------------

This toolbox solves smooth optimization problems defined on the manifold of Hierarchical Tucker (HT) tensors. There are also methods for interpolating HT tensors from missing entries. It is designed to solve the following problems

minimize_{ HT tensor X } 1/2 * ||P_{Omega} X - b||^2 + lambda R(X)

where P_{Omega} restricts to a specific set of indices of the tensor X, b is the subsampled data, lambda >= 0 and R(X) is a regularizer that prevents ranks from becoming degenerate. We refer to the References section for more details.

Curt Da Silva, March 2014

curtd@math.ubc.ca

-----------------------------------------------
## Features
-----------------------------------------------

- Efficient, dense linear algebra version for interpolating tensors
   - forms the full tensor at each iteration of the algorithm
   - serial and parallel versions
   - supports real + complex data

- Sparse linear algebra routines
   - implicitly forms the full tensor at each iteration - very minimal memory requirements
   - more efficient when the number of dimensions, individual dimension size is large compared to amount of data available
   - serial and parallel versions, embarassingly parallel + minimal communication
   - supports real data (can interpolate real + imaginary parts of data separately)

- Efficient methods for truncation of full data volumes to HT compressed form

- Tucker tensor-based interpolation methods, for comparison to Hierarchical Tucker methods

- Fully unit tested package 

-----------------------------------------------
## Installation
-----------------------------------------------

- For the current session - on the Matlab command line or in a script: 

   addpath(location/of/HTOpt/files);

- For every session - add the above line to your startup.m file

- To enable mex capabilities for tensor interpolation, run the installMex.m file in the main HTOpt directory

--------------------------------------------
## Examples
--------------------------------------------

Example interpolation scripts can be found in the /SLIM_ROOT/applications/processing/HierarchicalTuckerOptimization/examples directory. 

--------------------------------------------
## Directory structure
--------------------------------------------

examples/ - interpolation examples
- convergencespeed.m - plots the convergence speed of interpolating HT tensors for various optimization methods
- singlereflector_missingall.m - single reflector data interpolation for missing points
- singlereflector_missingrecs.m - single reflector data interpolation for missing receivers
- syntheticTest.m - generates interpolation experiments for synthetic data

htuck/ - dimension tree object + helper files
- dematricize.m - reshapes a matricized tensor back to a tensor
- dimensionTree.m - core dimension tree object
- dimTreeItr.m - abstract iterator object for iterating over a dimensionTree
   - dimTreeItrUp, dimTreeItrDown, dimTreeItrIntUp, dimTreeItrIntDown, dimTreeItrLeaves - various dimension tree iterators
- matricize.m - reshapes a tensor into a matrix
- truncate_ht.m - truncates a provided tensor in HT format to a specific accuracy + gives dimensionTree object
- ttm.m - matrix-tensor multiplication (multilinear product)
- ttt.m - tensor-tensor contraction

htuckOpt/ - optimization methods for the HT format
- gramian_regularizer.m - gramian-based regularizer for HT optimization (see paper for more details)
- minFunc_hTuck.m - main optimization routine (SD, CG, GN methods)
- opGramianJ.m - Gramian Jacobian SPOT operator
- opHTuckJ.m - HTuck derivative SPOT operator
- opHTuckJ2.m - same as opHTuckJ.m, except intermediate orthogonal projections omitted
- oppHTuckJ.m - parallel HTuck derivative SPOT operator
- oppHTuckJ2.m - parallel version of opHTuckJ2.m
- project_horizontal.m - projection on to the horizontal space at a point
- project.m - QR-orthogonalization retraction

interpolation/ - interpolation methods
- fitHT.m - wrapper for interpolating HT tensors
- HTInterpolation.m - convenience wrapper for subsampling + interpolating HT tensors
- ndimSubsampling.m - convenience function, makes n-dimensional subsampling operators

objective/ - misfit functions
- LSLinearMisfit.m - least-squares misfit for interpolation
- LSMisfitHT.m - sparse least-squares misfit for HT interpolation
- LSMisfitHTmex.m - sparse least-squares misfit for HT interpolation using mxlsmisfitht.c
- mxlsmisfitht.c - mex implementation of LSMisfitHT.m

spot/ - miscellaneous spot operators
- opPermute.m - permutation operator
- opProjection - subspace projection operator

tests/ - xunit tests
- testHTuck/ - unit tests for basic HT operations
- testHTuckOpt/ - unit tests for HT optimization operations
- testgeomCG/ - demonstrates equivalence of reference geomCG vs author-provided geomCG (Tucker tensor interpolation)
- testPackage.m - runs the tests from the above two directories

tuckOpt/ - tucker optimization methods (for comparison to HT methods)
- fitTucker.m - tucker interpolation
- fitTuckerRankIncrease.m - tucker interpolation w/ rank continuation
- HOSVD.m - higher order SVD

utility/ - miscellaneous functions
- cell_skeleton.m - creates a cell array copy of the input
- imagePlot.m - convenience function for displaying images
- idx1d2nd.m - convert 1d vectorized indices in to a matrix of Nd indices
- matrix2table.m - converts a matrix (of results) to a properly spaced table
- padarray.m - 2D/3D array padding function supporting singleton dimensions
- process_options.m - processes function options
- randsvd.m - randomized SVD for low-rank matrices
- qr_std.m - mathematically standard QR decomposition
- vec.m - vectorizes an ND array
- vec2TreeCell.m/treeCell2Vec.m - vectorizes/devectorizes quantities in cell array trees 


--------------------------------------------
## References
--------------------------------------------
Primary paper
C. Da Silva and F. J. Herrmann. "Optimization on the Hierarchical Tucker manifold - applications to tensor completion", SLIM Tech Report. Available at https://slim.gatech.edu/content/optimization-hierarchical-tucker-manifold-applications-tensor-completion 

Secondary papers
C. Da Silva and F. J. Herrmann. "Hierarchical Tucker Tensor Optimization - Applications to 4D Seismic Data Interpolation", EAGE 2013, available at https://slim.gatech.edu/content/hierarchical-tucker-tensor-optimization-applications-4d-seismic-data-interpolation-0

C. Da Silva and F. J. Herrmann. "Hierarchical Tucker Tensor Optimization - Applications to Tensor Completion", SAMPTA 2013, available at https://slim.gatech.edu/content/hierarchical-tucker-tensor-optimization-applications-tensor-completion