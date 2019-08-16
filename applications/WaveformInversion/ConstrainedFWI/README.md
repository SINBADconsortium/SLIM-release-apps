# Constrained FWI - regularizing nonlinear inverse problems using projections onto intersections of convex and nonconvex sets

##  DESCRIPTION
 This application is currently in development. Not all codes have input/output descriptions yet or propper comments. This application shows how constraints can be added to nonlinear inverse problem, while keeping codes to compute function values and gradients separated from optimization and constraints. The constraints and optimization are also separated. Although designed for nonlinear problems, this toolbox can also be used for linear problems.
 
 Suppose we want to minimize a data-misfit function ``f(m)`` which depends on medium parameters ``m: \min_{m} f(m)``. We define constraints as ``m \in C``, where ``C`` is usually a convex set. All algorithms describe here will also run if the set ``C`` is not convex, but it is not guaranteed to work well. ``C`` can describe bound constraints, various types of minimum smoothness constraints, convex and nonconvex transform domain sparsity constraints (TV, curvelet, etc.) among others. The set ``C`` is usually the intersection of multiple constraint sets. The user needs to provide a function handle `fh`, which takes the current model estimate as input and gives the function value ``f(m)`` and the corresponding gradient ``g`` as output. `[f,g]= fh(m)`. This package provides the tools to solve ``\min_{m} f(m) s.t. m \in C``.
 
##  COPYRIGHT
 You may use this code only under the conditions and terms of the
 license contained in the files LICENSE or COPYING provided with this
 source code. If you do not agree to these terms you may not use this
 software.
 
##  PREREQUISITES
 All prerequisites, except for MATLAB, are provided in the software
 release repositories and should be installed as necessary before using
 any of SINBAD's software.
 
##  INSTALLATION

###  Software in SLIM-release-comp repository
 In case this software repository is installed properly, running startup.m will load all the required paths. In case it does not run, add the following folders to the path:
 
 
 Folders which contain shared functions:

```code 
'/tools/solvers/SPGL1-SLIM'

'/applications/WaveformInversion/2DWRI/mbin'

'/tools/solvers/QuasiNewton/minConf_mod'

'/tools/solvers/SplittingMethods'

'/applications/WaveformInversion/2DWRI-TVconstrained'
```

The example script exp_2Dsimple_constrained.m also uses the following

```code
'/tools/algorithms/CommonFreqModeling'

'/tools/algorithms/2DFreqModeling'

'/tools/algorithms/3DFreqModeling'

'/tools/solvers/Krylov'
```

##  RUNNING
1: set up function handle `[f,g]= fh(m)`

2: Select which constraints should be used, for example

```matlab
	constraint.use_bounds       = 1;
	constraint.use_min_smooth   = 0;
	constraint.use_nuclear      = 0;
	constraint.use_rank         = 1;
	constraint.use_TV           = 1;
	constraint.use_cardinality  = 0; 
	constraint.use_TD_bounds    = 0; 
	constraint.use_cvxbody      = 0;
```

3: provide specifics for the selected constraints, for example:

```matlab
	%bounds
	constraint.v_max   = 5000; %upper bound
	constraint.v_min   = 1400; %lower bound

	%rank
	constraint.r_initial   = 2; %initial rank of the estimated model
	constraint.r_increment = 1; %increase the rank of the estimated model by this number each new frequency batch
```

4: obtain projectors onto each set `P{1}`, `P{2}`, ... each `P{i}` is a function handle which projects the input onto a set: `[projected_vector]=P{i}(input_vector)`. There are several scripts which set up the separate projectors (scripts with names like setup_constraints.m). These scripts should be modified and adapted for other applications.

```matlab
	P = setup_constraints(constraint,model_t,model_t,params,m0,1); %get projectors
```

5: get a function handle which projects the input onto the intersection of all selected constraint sets. By default this is (a modification of) Dykstra's algorithm.

```matlab
	funProj = @(input) Dykstra_prox_parallel(input,P,constraint.options_dyk); % get projector onto the intersection
```

6: use an optimization algorithm which handles projections properly. For example, the spectral projected gradient algorithm.
This example uses the minConf_SPG code from: [https://www.cs.ubc.ca/~schmidtm/Software/minConf.html](https://www.cs.ubc.ca/~schmidtm/Software/minConf.html) with licence [https://www.cs.ubc.ca/~schmidtm/Software/copyright.html](https://www.cs.ubc.ca/~schmidtm/Software/copyright.html)

```matlab
	[mn_t,fsave,funEvals,projects,iter_save] = minConf_SPG(fh,m0,funProj,options);
```

The example script implements this sequence of steps using the 2D frequency domain modelling algorithms from `SLIM-release-apps/tools/algorithms/2DFreqModeling/`

###  Preparing shell environment

 You must setup your shell environment according to the steps listed in
 the `README` located in home directory of the software release.

###  Downloading data

To run the example: 

```matlab
 run('SLIM-release-apps/applications/WaveformInversion/ConstrainedFWI/exp_2Dsimple_constrained.m'), 
``` 

no data needs to be downloaded, because it is a small example and data is generated on the fly.
 
###  Running applications/demos

```matlab
 run('SLIM-release-apps/applications/WaveformInversion/ConstrainedFWI/exp_2Dsimple_constrained.m')
``` 

requires parallel Matlab.

####  Hardware requirements

All codes which setup constraints, apply projections or do the optimization have serial versions. No parallel Matlab required for this. Some examples use our 2D frequency domain modeling which uses parallel Matlab.

##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
 
##  REFERENCES
This package includes the work described in:

1. [https://slim.gatech.edu/content/regularizing-waveform-inversion-projection-intersections-convex-sets](https://slim.gatech.edu/content/regularizing-waveform-inversion-projection-intersections-convex-sets)
2. [https://slim.gatech.edu/content/regularizing-waveform-inversion-projections-convex-sets-–-application-2d-chevron-2014-synthe](https://slim.gatech.edu/content/regularizing-waveform-inversion-projections-convex-sets-–-application-2d-chevron-2014-synthe)
