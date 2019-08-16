function w=FFTwaveletanalysis(x,FFTanalysisfilters,J);
% FFTWAVELETANALYSIS FFT-based implementation of the wavelet transform.
% 	w=FFTwaveletanalysis(x,FFTanalysisfilters,J) computes the wavelet 
% 	transform of a signal x using a Fourier method.  It uses periodic 
% 	boundary conditions.  The length of the signal must be a power of 
% 	two and the frequency responses of the filters specified in FFTanalysisfilters.
% 	
% 	Input:
% 	x=input signal, of size 2^N=length of FFTanalysisfilters
% 	J=depth of the decomposition, i.e., J wavelet bands (wav1 to wavJ) 
% 	+ 1 lowpass (lowJ) FFTanalysisfilters=[lowpassfilter;highpassfilter]
% 	
% 	Output:
% 	w=[wav1 wav2 ...  wavJ lowJ]: vector of pooled wavelet 
% 	coefficients, size 2^N
% 	
% 	See also FFTWAVELETSYNTHESIS, FFTFRACTSPLINEFILTERS, WEXTRACT.
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

M=length(x);
if M~=2^round(log(M)/log(2))
	disp(' ')
	disp('The size of the input signal must be a power of two!')
	disp(' ')
	w=[];
	return
end

if length(FFTanalysisfilters)~=M
	disp(' ')
	disp('The size of the input signal and of the filters must match!')
	disp(' ')
	w=[];
	return
end

% Fourier transform of the signal
X=fft(x,M);

G=FFTanalysisfilters(1,:);
H=FFTanalysisfilters(2,:);

w=[];
for j=1:J
	%
	% Computation of the outputs y and z
	%
	Y=G.*X;
	Z=H.*X;
	Y=1/2*(Y(1:M/2)+Y(M/2+(1:M/2)));
	Z=1/2*(Z(1:M/2)+Z(M/2+(1:M/2)));
	z=ifft(Z,M/2);
	w=[w z];
	
	M=M/2;
	X=Y;
	G=G(1:2:length(G));
	H=H(1:2:length(H));
end
w=real([w ifft(X,M)]);

