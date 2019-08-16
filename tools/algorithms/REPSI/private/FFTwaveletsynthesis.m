function x=FFTwaveletsynthesis(w,FFTsynthesisfilters,J);

% FFTWAVELETANALYSIS FFT-based implementation of the inverse wavelet 
% 	transform.
% 	x=FFTwaveletsynthesis(w,FFTsynthesisfilters,J) computes the inverse 
% 	wavelet transform of w.  This function is the inverse of 
% 	FFTwaveletanalysis.  It uses periodic boundary conditions.  The 
% 	wavelet coefficient vector has the following fine-to-coarse 
% 	organization: w=[wav1 wav2 ...  wavJ lowJ]
% 
% 	Input:
% 	w=wavelet transform, of size 2^N=length of FFTsynthesisfilters
% 	FFTsynthesisfilters=[lowpassfilter;highpassfilter]
% 	J=depth of the decomposition, i.e., J wavelet bands + 1 lowpass
% 	
% 	Output:
% 	x=signal of size 2^N
% 	
% 	See also FFTWAVELETANALYSIS, FFTFRACTSPLINEFILTERS, WEXTRACT.
% 	
% 	Author: Thierry Blu, October 1999
% 	Biomedical Imaging Group, EPFL, Lausanne, Switzerland.
% 	This software is downloadable at http://bigwww.epfl.ch/
% 	
% 	References:
% 	[1] M. Unser and T. Blu, "Fractional splines and wavelets," 
% 	SIAM Review, in press.
% 	[2] M. Unser and T. Blu, "Construction of fractional spline wavelet bases," 
% 	Proc. SPIE vol 3813, Wavelet Applications in Signal and Image 
% 	Processing VII, in press. 

M=length(w);
if M~=2^round(log(M)/log(2))
	disp(' ')
	disp('The size of the input signal must be a power of two!')
	disp(' ')
	x=[];
	return
end
if length(FFTsynthesisfilters)~=M
	disp(' ')
	disp('The size of the input signal and of the filters must match!')
	disp(' ')
	w=[];
	return
end

%
% Reconstruction of the signal from its
% bandpass components
%

G=conj(FFTsynthesisfilters(1,:));
H=conj(FFTsynthesisfilters(2,:));

M=M/2^J;
y=w(length(w)+((-M+1):0));
w=w(1:(length(w)-M)); 
Y=fft(y,M);
for j=J:-1:1
	z=w(length(w)+((-M+1):0));
	w=w(1:(length(w)-M));
	Z=fft(z,M);
	M=2*M;
	
	H1=H(1:2^(j-1):length(H));
	G1=G(1:2^(j-1):length(G));
	
	Y0=G1(1:M/2).*Y+H1(1:M/2).*Z;
	Y1=G1(M/2+(1:M/2)).*Y+H1(M/2+(1:M/2)).*Z;
	Y=[Y0 Y1];
end

x=real(ifft(Y,M));

