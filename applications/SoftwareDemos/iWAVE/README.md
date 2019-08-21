# Examples for using iWAVE interface for diffrent applications
##  DESCRIPTION
This application is an interface to use iWAVE++ in MATLAB. We provide easy functions to simulate the acoustic wave-equation, perform the born modeling and migration by using iWAVE++ from MATLAB. We provide simple scripts to use the forward modeling, linearized modeling and adjoint of linearized modeling. We also provide simple scripts to do RTM, least-squared Migration and FWI. This application only support the iWave in the Madagascar version 9520. In this new version, we provide functions for regularly sampled data, randomly sampled data and simultaneous sampled data. The usage of regular sampling does not require the support of Parallel matlab toolbox, while the usage of random sampling ang simultaneous sampling require the support of Parallel matlab toolbox. In this version, we only support the acquisition system that all shots are located at the same depth. We will release a version that supports source located at any arbitrary positon later. 

##  ON-LINE DOCUMENTATION
 <https://slim.gatech.edu/SoftwareDemos/misc/iWAVE/index.html>

##  COPYRIGHT
 You may use this code only under the conditions and terms of the
 license contained in the files LICENSE or COPYING provided with this
 source code. If you do not agree to these terms you may not use this
 software.

##  PREREQUISITES
Using this software requires installing iWAVE software (see INSTALLATION section below). The CWP Seismic Unix (http://www.cwp.mines.edu/cwpcodes/) is also required and it has to be compiled with XDR. All other prerequisites, except for MATLAB, are provided in the software
release repositories and should be installed as necessary before using
 any of SINBAD's software.

## INSTALLATION ##

###  Software in SLIM-release-comp repository
This application requires installation of extra software from SLIM-release-comp repository.

###  Software in SLIM-release-apps (this) repository
 Follow the instructions in the INSTALLATION file (located in the home
 directory of this software repository) to install necessary
 components.

###  Installing iWAVE
In order to use this application, you need to install iWAVE distributed with Madagascar version 9520 on your system. (You might want to check with your IT staff if the appropriate version of iWAVE is already installed on your system.) Here is the link that you can use to obtain it from SVN repository using the following command:
 
		svn co -r 9520 https://github.com/ahay/src/trunk/trip iWAVEv9520 
 
Installation of iWAVE depends on the type of MPI environment and compilers. Please, refer to instructions included in obtained source code.

#### Additional setup for iWAVE
The following elements must be added to shell's PATH environment:

		path_to_iWAVEv9520/rvl/seq/main
		path_to_iWAVEv9520/iwave/acd/main
		path_to_iWAVEv9520/iwave/asg/main
		path_to_iWAVEv9520/iwave/base/main
		path_to_iWAVEv9520/iwave/helm/main
		path_to_iWAVEv9520/iwave/trace/main
		path_to_iWAVEv9520/iwave++/acd++/main
		path_to_iWAVEv9520/iwave++/asg++/main

where *path_to_iWAVEv9520* above is an absolute path to home of your iWAVE installation.

###  Preparing shell environment
 You must setup your shell environment according to the steps listed in
 the README located in home directory of the software release.

##  DOCUMENTATION

You can use this application to simulate the time domain acoustic wave data by using iWave. You can also use this application to simulate the Born modeling data and its adjoint. You can also use this application to do time-domain RTM, least-square migration and full-waveform inversion. They are in separated demo folders ('ForwardModeling, LinearizedModeling, AdjointModeling, RTM, LSM, FWI').

The forward modeling, adjoint modeling and linearized modeling are three basic modeling operator. In order to use them, you need to provide a .txt file which contains necessary information for iWave correspondingly and give the file name to 'options.fwdpara, options.adjpara and options.linpara'. In the demo, we only provide some basic parameters that you may need. If you want to know all parameters, please read files in the directory './iWaveParameters'.

We provide script to perform the forward modeling, linearized modeling and adjoint modeling with regularly sampled data, randomly sampled data and simultaneously sampled data. You should be able to find scripts in the directory './demo/ForwardModeling', './demo/LinearizedModeling', and './demo/AdjointModeling' with the following script files:

1. 'Forward_Demo.m' - Forward modeling script with regularly sampled data;

2. 'Forward_Demo_Rand.m' - Forward modeling script with regularly sampled data and simultaneously sampled data;

3. 'Linear_Demo.m' - Linearized modeling script with regularly sampled data;

4. 'Linear_Demo_Rand.m' - Linearized modeling script with regularly sampled data and simultaneously sampled data;

5. 'Adjoint_Demo.m' - Adjoint of linearized modeling script with regularly sampled data;

6. 'Adjoint_Demo_Rand.m' - Adjoint of linearized modeling script with regularly sampled data and simultaneously sampled data;

When you use regular sampling, **you do not need to open matlabpool**. The variable `options.nlab` defines the total number of processes you want to use. It equals to the product of `partask`, `mpi_np1` and `mpi_np2` in the parameter `.txt` file.

When you use random sampling and simultaneous sampling, **you should open matlabpool, and the number of matlabpool you open equals how many shots you want to simulate at the same time**. In this case, the variable `options.nlab` defines the number of processes you want to use for simulating each shot. It equals to the product of `mpi_np1` and `mpi_np2` in the parameter `.txt` file. And the variable `partask` in the parameter `.txt` file should be `1` in this case. If you are using multiple nodes, please distribute the matlab workers equally into these nodes. And make sure the product of number of workers in each node and `options.nlab` is not larger than the total number of processes in each node.

In the scripts with randomly sampled data and simultaneously sampled data, you may need to modify the parameter 'flagrandseq' and 'flagsim' to decide which sample scheme you want to use. 'flagrandseq=1' means that you select random sampling. 'flagsim=1' means that you select simultaneous sampling. You should not set these two parameters be 1 at the same time. For simultaneous sampling, you can set the parameter 'nsupshot' which defines how many simultaneous shot you want to generate.  


##  RUNNING
 Before you run any one of scripts, you need to know the mpirun command in your environment and change it in all scripts. You need to change the content inside the file : `tools/utilities/iWAVE/iwave_MPIcommand.m`. Replace the sentence `str = ['mpiexec.hydra -np ' num2str(nlab) ];` to your own mpirun command in your environment, for the usage of regular sampling. And also Replace the sentence `str = ['mpiexec.hydra -np ' num2str(nlab) ' -hosts localhost'];` to your own mpirun command with local hostfile in your environment for the usage of random sampling and simultaneous sampling. Two examples are in the file. You need to run the 'startup.m' in this application first to add all necessary toolbox. 
 
 After changing the file `tools/utilities/iWAVE/iwave_MPIcommand.m`, you can directly run all scripts in the './demo' directory without changing anything.

###  Preparing shell environment
 You must setup your shell environment according to the steps listed in
 the README located in home directory of the software release.

###  Downloading data
 There are no data required to be downloaded for this application. But you need to run the script 'Generate_Model_data.m' in the director './demo_model' to generate all necessary data and model for demo scripts.

####  Hardware requirements
 Number of CPU cores required: 4 or 16
 
 Amount of memory required: 4 GB
 
 Expected runtime: 
 
 for applications except FWI,LSM and RTM, less than 2 minutes
 
 for applications LSM (3 iterations) and RTM, less than 20 minutes
 for applications FWI, about 1 hour
 

####  Data adaptation
If you want change the number of processes, you need to change both the demo file and parameter file. In the demo file, you need to change the value of 'options.nlab' to the number of processes you want to use. In the parameter.txt file, you need to change the value partask (number of shots simulated parallelly), mpi_np1 (number of parts that the first direction will be devided) and mpi_np2 (number of parts that the second direction will be devided). ``n_{process} = partast * mpi_np1 * mpi_np2``. 

If you want to use your own data, please make sure that the order of your data should be 'time / receiver / source'.

The model parameter used in this application is bulk modulus (= velocity(km/s)^2*density(``g/cm^3``)) and buoyancy (1/density(``g/cm^3``)).              

##  NOTES
 Please make sure that the demo directory is not too deep. Otherwise, you may need to copy the demo directory to some shallow place such as 'Desktop' and change all paths in demo scripts for necessary toolbox, data and model.
 
 Right now, the gradient of FWI has some very big value at the source location, so if you want to use this application to do FWI for your own job, you'd better make a mask to mute the gradient around the source locations. The 'misfit_cut.m' in the FWI is an example.
 
 This version, we only support the case that all shots and receivers are alined with constant interval. Randomly selected shot and simultaneous shot are not supported.
 
##  SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.

##  REFERENCES
 1. Symes, W.W., Sun, D. and Enriquez, M., 2011. From modelling to inversion: designing a well‐adapted simulator. Geophysical Prospecting, 59(5), pp.814-833.
 2. Tschannen,V., Fang, Z., and Herrmann, F. J., 2014. Time domain least squares migration and dimensionality reduction, Technical report, UBC.
