FRACTIONAL SPLINE WAVELET SOFTWARE FOR MATLAB

Thierry Blu and Michael Unser
Biomedical Imaging group 
Swiss Federal Institute of Technology Lausanne (EPFL), Switzerland.
October 1999

A. TEST PROGRAM

DEMOFRACT: Demonstration program that generates fBm-like noise using 
orthogonal symmetrical fractional spline wavelets.  For example, try J=12, 
alpha=0.3, beta=1.3.  You can modify this file to use non-symmetrical 
fractional splines, or also semi-orthonormal fractional splines (modify the 
'type' variable: type = '+ortho', or '-ortho', or '*ortho', or '+bspline', or 
'-bspline', or '*bspline', or '-dual', or '+dual', or '*dual')

DEMOFRACT2: Generates a series of fBm-like noises.


B. BASIC WAVELET TRANSFORM FUNCTIONS

FFTFRACTSPLINEFILTERS: Computes the frequency responses of the wavelet 
analysis and synthesis filters (lowpass and highpass) for various types 
('+ortho','-ortho','*ortho','+bspline','-bspline','*bspline','-dual','+dual','*dual') 
of orthonormal or semi-orthonormal fractional spline wavelet transforms.  
The type "*ortho" is a continuous-order generalization of the 
Battle-Lemarié transform.

FFTWAVELETANALYSIS: Performs a FFT-based computation of the wavelet 
transform of a signal to a given depth J.

FFTWAVELETSYNTHESIS: Performs a FFT-based computation of the inverse 
wavelet transform.


C. OTHER ROUTINES

WEXTRACT: Extraction of a component (or band) of the wavelet 
transform.

FRACTSPLINEAUTOCORR: Frequency domain computation of fractional spline 
autocorrelation function.





