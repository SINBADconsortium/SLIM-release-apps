#3D full-waveform inversion
This applications is available only in the software release for members of SINBAD consortium.

#DESCRIPTION
This package contains a matlab scripts that illustrate the 3D waveform
    inversion algorithm as described in the Technical Report 
    "3D Frequency-domain seismic inversion using a row-projected Helmholtz solver",
    T. van Leeuwen and F.J. Herrmann, 2013.<br />
    An overview of the results can be found in the ./doc directory and online:<br />
    <https://slim.gatech.edu/SoftwareDemos/applications/WaveformInversion/3DBatching/>
#COPYRIGHT
You may use this code only under the conditions and terms of the
    license contained in the file LICENSE or COPYING.txt provided with
    this source code. If you do not agree to these terms you may not
    use this software.
#PREREQUISITES
All prerequisites, except for MATLAB, are provided with the
    software release and should be installed before using any of
    SINBAD's software.
#INSTALLATION
Follow the instructions in the INSTALLATION file (located in the
    root directory of this software release) to install all 3-rd party
    software (except for MATLAB) and SINBAD's software.
#DOCUMENTATION
Documentation and examples in .html format are included in the ./doc
subdirectory.  Help can also be accessed from Matlab via the `help' command.

#RUNNING
Start matlab from this directory or run startup.m to set the correct paths.

1. To run the scripts, first download the data and models from<br />
        <ftp://ftp.slim.gatech.edu/data/SoftwareRelease/WaveformInversion/3DBatching/data/>
       by typing `scons' in the ./data directory
2. The scripts that produced the examples shown in the documentation can be found in ./scripts.

#NOTES

## Controlling the number of threads

Some linear solvers in this application (such as `(p)CARP(B)CG.m`) make use of
a MEX file: `sweepR_mex`. This MEX file implements the Kaczmarz (serial) or
CARP (parallel) algorithms. Parallelism is implemented via threads
(specifically, the `pthreads` library). By default, the number of threads is
set to one. It can be controlled in two ways: either via the `OMP_NUM_THREADS`
environment variable or via an optional extra argument to `sweepR_mex`. To
specify the number of threads via an environment variable, execute the
following at the MATLAB prompt, replacing `'10'` with the desired number of
threads:

```
setenv('OMP_NUM_THREADS', '10');
```

To use instead the extra optional argument to `sweepR_mex`, call it as follows:

```
sweepR_mex(R,idx,x,b,w,d,10);
```

Where the first 6 arguments are the mandatory ones, explained in the header
comment of the source code file `sweepR_mex.c`, and the last argument is the
number of threads to use. 

NOTE: the optional argument, if provided, takes precedence over the environment
variable.  

If there are multiple MATLAB workers on one node, be careful to not set the
OMP_NUM_THREADS too high, since each worker doing parallel CARP sweeps will use
that many threads. For example, setting OMP_NUM_THREADS to 32 while there are 4
MATLAB workers per node means there will be 32 x 4 = 128 parallel CARP threads,
which may be too high. Experiment to find the best combination of MATLAB MPI
parallelism and parallel CARP sweeps via multi-threading on your machine.

#SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
