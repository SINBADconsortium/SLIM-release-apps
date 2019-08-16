#Sequential interpolation using weighted L1 minimization
##DESCRIPTION
This package contains Matlab functions for interpolation of randomly subsampled data volumes by solving weighted L1 minimization sequentially on partitions of the data volume.

    b         - subsampled data volume
    mask      - logical data volume indicating the locations of the missing data entries
    options:
	partOrder - specifies the proper permutation of the data volume so that partitioning is always performed along the first dimension
	st        - index of the first partition
	fin       - index of the last partition
	C	  - sparsifying transform
	maxiter   - maximum number of spgl1 iterations
	omega	  - sets the weight used in weighted L1 minimization
##COPYRIGHT
You may use this code only under the conditions and terms of the
    license contained in the file LICENSE or COPYING.txt provided with
    this source code. If you do not agree to these terms you may not
    use this software.
##PREREQUISITES
All prerequisites, except for MATLAB, are provided with the
    software release and should be installed before using any of
    SINBAD's software.
#INSTALLATION
Follow the instructions in the INSTALLATION file (located in the
    root directory of this software release) to install all 3-rd party
    software (except for MATLAB) and SINBAD's software.
##DOCUMENTATION
Documentation for each of the functions can be accessed by typing `help <function>' in Matlab. 
    Examples are provided in applications/Interpolation/WeightedL1 
##RUNNING
The functions can be called directly from Matlab. 
##NOTES
##SUPPORT
 You may contact SLIM's developers of SINBAD software via issue tracker for this repository. We do not have resources to actively support public version of our software. However, we will try to answer the questions as much as possible.
