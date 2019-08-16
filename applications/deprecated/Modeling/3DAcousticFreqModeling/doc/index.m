%% 3D Frequency-Domain Modeling Kernel
%%  Description
%
% This package contains the Matlab functions to manipulate the discrete 3D 
% acoustic wave equation operator in frequency domain. It is primarily based on
% a staggered grid dispersion minimizer second order finite differences stencil
% first proposed in [1], using PML. In addition, we currently have implemented 
% the standard second order 7 points finite differences stencil with PML 
% (c.f. [2] for more details).
%
%% Running & Parallelism
% All the examples can be reproduced by running the scripts found in
% the software release under |applications/Modeling/3DAcousticFreqModeling|. 
% Start matlab from that directory or run |startup| in that directory to add 
% the appropriate paths.
%
% To properly run some of the demos and tests, it is necessary to download data
% files from the SLIM FTP server. 
%  
% *Installation with SCONS*
% 
% Enter the data directory and use the scons executable:
%      
% |cd ./data|
% 
% |scons|
%      
% *Installation without SCONS*
% 
% Download all the data files found in 
%  
% |ftp://ftp.slim.gatech.edu/data/SoftwareRelease/Modeling/3DAcousticFreqModeling/data|
%      
% and place them into ./data folder. The files are |m_init.rsf|, |m_true.rsf| 
% and |Dobs.rsf|.
% 
% 
% *Remark:* Currently, most of this code is for sequential use only. Most 
% notabily, it does not discretize the Helmholtz operator in parallel, but in 
% sequential mode. It is possible to distribute the discrete operator to then
% be used by parallel solvers such as CARPCR.
% 
% In future releases, this issue will be tackled.
%
%
%% Functions 
% 
% Follows a detailed description of all the functions provided in this toolbox.
% Use 
% 
% |help function_name|
%
% for more details on the input and output of each function. Unless otherwise 
% specified, all function in this package are in its stable release, meaning 
% that its efficiency in computational time, accuracy and memory use has been 
% confirmed through several numerical experiments. Those marked as (beta) and 
% (alpha) should function properly, but have not been extensively tested yet and
% are thus subject to improvements in the future.
%  
% * *WaveEquationFD (alpha)*
% 
% The main function of this package. Given a velocity model and an acquisition
% geometry with its respective dataset, this returns the predicted data, the 
% value of the misfit function for that dataset and the gradient. Each output
% parameter is optional and can be disabled.
% 
% This function handles the whole discretization process; it shrinks/stretches
% the domain to satisfy stability conditions depending on the frequency. It also
% adds the PML of the proper size and gives an option to cut the PML out or 
% to stack PML points (its adjoint). The full description of all the options
% for this function can be found in its own documentation. Use |help WaveEquationFD|
% for more details.
% 
% Example: <demo_WaveEquationFD.html examples/demo_WaveEquationFD.m>.  
%
%
% * *ForwardFD (beta)*
% 
% A wrapper to WaveEquationFD which returns only the computed data.
% This function (that is, the data output) has been tested much more extensively
% than WaveEquationFD and thus is considered almost in its stable release.
% Consult the demo for a comparison between the old and new forward modelling
% kernel using this function.
% 
% Example: <demo_ForwardFD.html examples/demo_ForwardFD.m>.
% 
% * *MisfitFD (alpha)*
% 
% A wrapper to WaveEquationFD which returns only the most vanilla misfit value
% (no regularization are used in this version) as well as the value of the
% gradient. This is the least tested component of WaveEquationFD function.
% NOTE: serial only at the moment.
% 
% Example: <demo_MisfitFD.html examples/demo_MisfitFD.m>.  
% 
% * *helmholtz_3d*
%
% Computes the discrete Helmholtz operator using a staggered grid second order
% finite differences scheme. The output is an array with  dimension 
% (27,nx*ny*nz), or in other words, it stores the discrete operator in a band-
% storage format rather than as a sparse matrix. Each column of this array 
% represents a point in the physical domain.
% This stencil is discussed in [1] and it aims at minimizing the disperson 
% error, attaining stability with as little as 4 points per wavelength in 
% homogeneous testcases. We usually recommend the use of 6 points per 
% wavelength for heterogeneous media.
% 
% Example: <demo_helmholtz_3d.html examples/demo_helmholtz_3d.m>.  
% See also <demo_R_Helm3D.html examples/R_Helm3D.m> for an example using the legacy code.
%
% * *discrete_helmholtz* 
% 
% Given a structure model (containing the velocity model and grid spacing) and
% the frequency, this function returns a structure containing all the parameters
% needed to compute the discrete operator with a stable computational grid 
% spacing. It automatically sets the PML size, interpolates the given model
% to the required size, and returns function handlers that allow the user to
% easily manage data both on the model grid and on the computational grid. It
% also returns the discrete operator (by default, but can be disabled to merely
% compute the proper stable parameters).
% 
% To understand better the output of this function, it is important to 
% distinguish the computational grid from the physical grid.
% 
% The physical grid is the grid in which the velocity model is provided. For 
% computing a numerical solution for the wave equation, however, we might need
% a finer (or coarser) grid. This is associated with the wavelength of the waves
% being propagated in the frequency domain, and thus depends on the target 
% frequency. Also, the computational grid may contain an artificial absorbing
% layer to avoid reflections.
% 
% Example: <demo_discrete_helmholtz.html script/demo_discrete_helmholtz.m>
%
% * *helmholtz_3d_derivative*
%
% Computes the tensor ∂H/∂m, where H is the Helmholtz matrix discretized using 
% helmholtz_3d.m routine, and m is the model. This is required, for instance,
% to compute the gradient of the misfit function efficiently in the adjoint-
% state formulation, instead of computing the full Jacobian. Although this is a
% tensor, it can be stored exactly as H itself, in band-storage format. Using
% this object might be a bit trickier than expected. Please, refer to the full
% documentation, the demos or Matlab's Help on how to use this tensor. Also
% a brief explanation on how to deduce this implementation can be found in [3].
%  
% * *Hmvp*
% 
% Multiplies a vector by a matrix in band-storage format (such as the operator
% generated with helmholtz_3d and helmholtz_3d_7p). It calls a mex  function to
% perform the matrix-vector product itself in an efficient and multi-threaded 
% environment.
%
% * *create_pml*
% 
% Function used exclusively to compute the damping coefficients for the PML 
% layer. This function is called internally and never needs to be explicitly
% invoked by the user.
%
% * *H2sparse*
% 
% Converts the output of helmholtz_3d to a standard Matlab's sparse matrix 
% format. Specially useful if you intent to use direct solvers.
%
% * *Htransp*
% 
% Returns the transpose of the input matrix, given in band-storage format.
% Notice that this does not conjugate the operator!
% 
% * *helmholtz_3d_7p (beta)*
% 
% Very similar to helmholtz_3d.m, except that this generates a discrete 
% operator using the standard second order 7 points stencil with the PML; 
% Please refer to [2] Appendix A for more details on this stencil. 
% This generates a  matrix with 7 diagonals only (as opposed to 27) but 
% requires more points per wavelength. In the literature it is common the use 
% of 10 points per wavelength, but 12 points is a more conservative and also 
% a popular choice.
%
% * *helmholtz_solver (pre_alpha)*
% 
% Replaces the old setup_solver and partially replaces discrete_helmholtz
% function. It follows the same syntax as discrete_helmholtz but returns the
% discrete operator as wel as a pointer for the linear solver. The solver is
% chosen based on the size of the problem, discretization used, amongst other
% characteristics of the problem. The user can specify target tolerance, 
% maximum number of iterations, or even to impose a specific solver. Please,
% refer to "help helmholtz_solver" for  more details about the interface.
% 
% IMPORTANT: the use of direct solvers is relatively crippled (by the fact 
% that it does not store the incomplete factors, but re-compute them for every 
% single right hand side). This function should see large improvements in the 
% forthcoming releases. 
%
% *Remarks for Developers:*
%
% # This package allows free surface; this is achieved by setting the PML length 
%   in the top of the domain to zero. See the examples for more details.
% # The legacy routines generates a (n_x*n_y*n_z, 27); this code inverts this to 
%   the more efficient (27, n_x*n_y*n_z).
% # The legacy version uses slowness squared. This code allows the user to 
%   choose between meters per second or slowness squared (this is emphasized in 
%   the  comments inside the code itself).
% # The legacy code permutes the z direction to the first dimension. This code 
%   keeps z direction as the third direction.
%

%% Tests
%
% *<test_validate.html test_validate.m>*
% 
% This test certifies that the discrete operator computes a satisfactory 
% numerical solution with acceptable dispersion error. This test suit also 
% compares the old discretization code with the new one. The Matlab file can be 
% found in the scripts directory.
%
% *<test_helmholt_derivative.html test_helmholt_derivative.m>*
% 
% This test guarantees that the tensor (stored in a band-storage format) 
% computed by helmholt_3d_derivative is in fact the derivative of the matrix
% generated by helmholt_3d.
%
% *<test_gradient.html test_gradient.m>*
% 
% To guarantee that the gradient computed by MisfitFD is in fact the numerical
% gradient of the misfit function, we written this script showing that as we 
% decrease the magnitude of the model perturbation δm, the gradient behaves as
% a "second order correction". 
% 


%% References
% <http://dx.doi.org/10.1190/1.2759835 [1]> Operto, S., J. Virieux, P. Amestoy, 
% J.-Y. L’Excellent, L. Giraud, and H. B. H. Ali, 2007, 3D finite-difference 
% frequency-domain modeling of visco-acoustic wave propagation using a 
% massively parallel direct solver: A feasibility study: Geophysics, 72, no. 5, 
% SM195–SM211
% 
% <http://ethesis.inp-toulouse.fr/archive/00001221/01/pinel.pdf [2]> Pinel, X.,
% 2010, A perturbed two-level preconditioner for the solution of 
% three- dimensional heterogeneous Helmholtz problems with applications to 
% geophysics, INPT PhD Thesis
% 
% <2014-Lago-Efficient_computation.pdf [3]> Lago, R. 2014 - Efficient 
% computation of dH/dm tensor for a 27-points stencil with mass lumping
%
