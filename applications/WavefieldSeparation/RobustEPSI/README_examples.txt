# Author      : Tim Lin
#               Seismic Laboratory for Imaging and Modelling
#               Department of Earth & Ocean Sciences
#               The University of British Columbia
#         
# Date        : Faburary, 2011
# 
# You may use this code only under the conditions and terms of the
# license contained in the file LICENSE provided with this source
# code. If you do not agree to these terms you may not use this
# software.

Robust Estimation of Primaries with Sparse Inversion (via L1 minimization)
         for Pre-Stack Seismic Lines
==============================================================================

    The main driver of Robust EPSI is 'EPSI_SLIM_main.m' nested in the 'tools'
    directory. The different folders here contain different processing
    examples/experiments of using the main driver. They each contain a SCONS
    makefile that produces individual Matlab scripts for different
    settings/dataset.

FETCHING DEMO SEISMIC DATA

    Before the demo could be run, go to the 'data' directory and run 'scons' 
    to retrieve (with cURL) the three demo seismic datafiles used in this 
    release. 
    
    If for some reason the fetching script failed to work, you can download a 
    zip archive of the files at the address below and unpack it into the 
    'data' folder:
    
    ftp://ftp.slim.gatech.edu/data/users/timlin/RobustEPSI/RobustEPSI_dataset.zip
    
    
RUNNING THE EXAMPLES
	
	After the datafiles have been retrieved:
	
	1)  - make sure SCons >=1.2 is installed
	    - return to the 'Experiments' directory.
	    - cd to the 'Synthetic' directory.
	    - run 'scons'
	
	2)  The 'Synthetic' directory contains settings to process a synthetic 
	    marine dataset with a salt-dome structure. Running 'scons' in here
	    generates three Matlab scripts:
	        
	        'Synthetic.m'
	        Run this script in Matlab to see Estimation of Primaries with L1 
	        Inversion without using sparsity-promoting transforms to help in 
	        estimating the Green's functions. (approx. runtime 1.5 hours, 2GB 
	        memory)
	        
	        'Synthetic_sparsity.m'
	        Same as above but with sparsity-promoting transform enabled. 
	        (approx. runtime 3 hours, 8GB memory)
	        
	        'Synthetic_sparsity_parallel.m'
	        Same as above but carries out the computation on remote computers 
	        using the Parallel Computing Toolbox in Matlab. Must be in a valid 
	        parallel pool session before running.
	        
	    For all these scripts, the solution (primary wavefield and estimated 
	    wavelet) will be written to a .mat datafile in the same directory, 
	    with obvious names.
	    
    3)  The SConstruct file is interpreted by Python (vis SCons) before 
        building the Matlab scripts. It is therefore possible (encouraged   
        even) to edit parameters directly inside the script.
    
        As an example, try changing the 3 experiments in the script to output 
        BIN files instead of .su files. This is done by simply changing the 
        output filename's extension. The software will automatically determine 
        output data format depending on the extension. Execute 'scons' again 
        to re-generate the new Matlab scripts and re-run them. Other changes
        you can try is setting 'verbosity' to 1, to see internal per-iteration
        output of the underlying solvers used.
        
        Notice how we are using native Python structures (dictionary, which is
        Python's implementation of a hashtable/struct) to define the
        experiments. Here the native Python command 'copy' is used to clone a
        common options set for all three versions of the processing. You can
        take advantage of this in may ways, for example, you can quickly use a
        for loop to implement a quick parameter search over a particular
        parameter.
        
    4)  The 'GulfOfSuez178' directory contains settings to process a real, 
        preprocessed shallow-ocean bottom marine dataset from the Gulf of 
        Suez. The setup is exactly the same as the 'Synthetic' folder. In this 
        case the datafiles are in the SU format.
    
        The 'GulfOfSuez355' directory contains the same experiment as above
        but uses a more finely sampled dataset.


USING YOUR OWN DATAFILES
        
        REPSI currently deals with 3 file-types: MAT, SU, and Binary Cubes
        (BIN). It uses the filename extension (.mat, .su, .bin respectively)
        to distinguish between these types. Binary cubes can be obtain, for
        example, by running 'sustrip' on a SU file, or from the datafiles of
        any dictionary-based filetypes.
        
         (NOTE: The modified SegyMAT library used for SU file i/o currently
        only handles native files. XDR support is currently not available. If
        you installation of SU is compiled with XDR support, please convert
        the datafiles to native format before using with REPSI. Similarly, you
        should also run commands like 'suoldtonew' to convert Matlab outputs
        back to XDR format.)
        
         The core computation involved in REPSI is based on the same
        multidimensional data-convolution that SRME performs for multiple
        prediction. It therefore theoretically also requires a split-spread,
        regularly sampled geometry. The implementation of REPSI is very
        simple, and therefore preprocessing to get data into the format of
        split-spread, regularly sampled shot gathers, with source and
        receivers lying on the same grid, is required. For now, missing
        near-offsets should also be interpolated, until simultaneous data
        interpolation is implemented in the current code.
        
         As such, note that special formatting rules apply to the input
        datafiles. The fast dimension (trace dim) is time, the second
        dimension is trace number, and the third dimension is the shot number.
        Traces must be constant length. Receiver position and shot position
        must range over the same values, implying that the second and third
        dimensions must have the same size. On some datasets (i.e., marine
        data) this requires some redatuming and/or interpolation. When taking
        a constant time slice, the image should be (nearly) diagonally
        symmetric due to reciprocity.
        
         If at all possible, decompose the data into up and down-going
        wavefields, and use the up-going wavefield with REPSI. However, this
        is often not critical. If the data has been de-ghosted, you should try
        setting useOblique to 1
        
         The option 'test_readData' can be specified to the main driver
        function to quickly check whether the dataset is successfully read by
        Matlab. Specify the 'input_endian' option to deal with possible endian
        issues.
        
        

SUMMARY OF APPROXIMATE SYSTEM REQUIREMENTS AND RUNTIME

    MATLAB currently requires all data to be stored in the system memory. As
    this demo is written in Matlab, the following resource requirements apply
    when running the demos.
    
    (NOTE: these are approximate process memory requirements that do not
    include resident system memory usage, and are only supplied as a
    guideline.)
    
    Synthetic/
        'Synthetic.m'
                - 1.5 Hours, 4 GB RAM
        'Synthetic_sparsity.m'
                - 3 Hours, 16 GB RAM
        'Synthetic_sparsity_parallel.m'
                - runtime varies based on available computation workers
        
    GUlfOfSuez178/
        'GulfOfSuez178.m'
                - 4 Hours, 7.5 GB RAM
        'GulfOfSuez178_sparsity.m'
                - 15 Hours, 24 GB RAM
        'GulfOfSuez178_parallel.m'
                - runtime varies based on available computation workers, total
                  8GB RAM
        'GulfOfSuez178_sparsity_parallel.m'
                - runtime varies based on available computation workers
    
    GUlfOfSuez355/
        'GulfOfSuez355.m'
                - 8 Hours, 24 GB RAM
        'GulfOfSuez355_sparsity.m'
                - 30 Hours, 120 GB RAM
        'GulfOfSuez355_parallel.m'
                - runtime varies based on available computation workers, total
                  24GB RAM
        'GulfOfSuez355_sparsity_parallel.m'
                - runtime varies based on available computation workers
                
