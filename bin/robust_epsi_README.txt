# Author: Tim T.Y. Lin
#         Seismic Laboratory for Imaging and Modeling
#         Department of Earth & Ocean Sciences
#         The University of British Columbia
#         
# Date  : Feb, 2011

Robust Estimation of Primaries with Sparse Inversion (via L1 minimization)
===================================================================================================

Documentation concerning the 3 executable files:

        robust_epsi
        robust_epsi_spars
        robust_epsi_parallel
        robust_epsi_parallel_spars

The four programs here are each Python wrappers to the main Matlab
driver 'EPSI_SLIM_main.m'. They REQUIRE that Matlab be able to be evoked
from your shell, in command-line only mode, by running:

    'matlab -nodesktop -nosplash'

These scripts can help you to run Stabilized EPSI on your own datasets, with your
own custom parameter settings:
    
    'robust_epsi'
        This is the vanilla Estimation of Primaries with L1 Inversion.
        
    'robust_epsi_spars'
        This version enables the use of sparsity-promoting transforms (wavelet 
        and curvelets). The ns_per_trace of the input data must be powers of 2.
        
    'robust_epsi_parallel'
        This version should only be used in conjunction with the Parallel 
        Computing Toolbox (v4.2+) on Matlab. A default configuration name must 
        be set by hand inside the script.
    
    'robust_epsi_parallel_spars'
        This version enables the features from the previous 2 programs.
        
Please look at the source codes to see the list of available and mandatory
options when invoking these programs.

The input and output files may have extensions '.su' (SU files), '.mat'
(Matlab datafiles), or '.bin' (float32 binary blobs at the machine's native
endian). The programs will automatically handle the datafiles correctly based
on their extensions, and check for possible endian issues. The '.bin' format
in particular is meant to be used with 'sustrip' and 'supaste' to workaround
issues with bad trace headers.

(NOTE: The modified SegyMAT library used for SU file i/o currently only handles
native files. XDR support is currently not available. If you installation of SU
is compiled with XDR support, please convert the datafiles to native format
before using with REPSI. Similarly, you should also run commands like 
'suoldtonew' to convert Matlab outputs back to XDR format.)

NOTE ON DATAFILES:
The core computation involved in REPSI is based on the same multidimensional 
data-convolution that SRME performs for multiple prediction. It therefore
theoretically also requires a split-spread, regularly sampled geometry. The
implementation of REPSI is very simple, and therefore preprocessing to get
data into the format of split-spread, regularly sampled shot gathers, with
source and receivers lying on the same grid, is required. It CANNOT, for example,
work directly on marine geometry. For now, missing near-offsets should also be 
interpolated, until simultaneous data interpolation is implemented in the 
current code.

As such, note that special formatting rules apply to the input datafiles. The
fast dimension (trace dim) is time, the second dimension is trace number, and
the third dimension is the shot number. Traces must be constant length.
Reciever position and shot position must range over the same values, implying
that the second and third dimensions must have the same size. On some datasets
(i.e., marine data) this requires some redatuming and/or interpolation. When 
taking a constant time slice, the image should be (nearly) diagonally symmetric
due to reciprocity.

EXAMPLE USAGE

    'robust_epsi -h'
    
        Displays help string and all available options
        
    'robust_epsi --input_file=input.su --test_dataread'

        Launches a test to see whether Matlab is able to read the input
        datafile 'input.su' in your current directory without issue.
        
    'robust_epsi --input_file=input.su --output_primary_file=output.su 
        --output_wavelet_file=output_wavelet.su --topmuteT=0.1'
        
        Runs the standard Stabilized EPSI (without sparsifying transforms) on
        an input SU file called 'input.su', writing out the primary wavefield
        into a file called 'output.su' and the estimated wavelet into a file
        called 'output_wavelet.su'. A muting mask of 0.1s from t=0 is applied
        to each trace on the predicted primary to prevent trivial solutions.
        
        All these options here are mandatory, except for the preceding two
        usage scenarios.
        
    'robust_epsi --input_file=input.su --output_primary_file=output.su 
        --output_wavelet_file=output_wavelet.su --topmuteT=0.1
        --input_endian=b'

        Same as above with the input datafile specified to be in big endian.
        The input_endian flag is used to solve potential endian issues.
    
    'robust_epsi --input_file=input.mat --output_primary_file=output.mat 
        --output_wavelet_file=output_wavelet.mat --topmuteT=0.1
        --q_estlength_posT=0.2 --q_estlength_negT=0'
        
        Same as the above, except the input and output files are now presumed
        to be in native Matlab form (input wavefield datacube must be in a
        variable called 'data', and a number 'dt' in seconds). 
        
         The estimated wavelet is inforced to be entirely
        causal (s_estlength_negT=0) with an allowed width of 200ms.
        Matlab native files do not have endian issues so the input_endian
        flag should not be needed.
        
    'robust_epsi --input_file=input.mat --output_primary_file=output.mat 
        --output_wavelet_file=output_wavelet.mat --topmuteT=0.1
        --window_startT=0.1 --window_endT=0.3'
        
        Run Stabilized EPSI but forces the first primary events to be picked
        from a time window that starts at 0.1s and ends at 0.3s for all
        traces. May be useful if the dataset is extremely noisy or if there
        are significant residual energy before the first reflected event.
        (Although topmuteT may also be used for this purpose.)
