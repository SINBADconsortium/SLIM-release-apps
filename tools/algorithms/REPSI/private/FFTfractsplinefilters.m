function [FFTanalysisfilters,FFTsynthesisfilters]=FFTfractsplinefilters(M,alpha,type)
% FFTFRACSPLINEFILTERS Filters for the various orthogonal fractional 
% 	spline wavelet transforms (orthonormal or semi-orthonormal).  
% 	[FFTanalysisfilters,FFTsynthesisfilters]=FFTfractsplinefilters(M,alpha,type) 
% 	Computation of the frequency response of the low- and high-pass 
% 	filters that generate the orthonormal or semi-orthonormal (B-spline  
%	or dual) fractional splines of degree alpha, and of type '+' (causal),
% 	 '-' (anticausal) or '*' (symmetric).
% 	
% 	Input:
% 	M : size of the input signal = length of FFTfilters = 2^N
% 	alpha>-0.5 : degree of the fractional splines 
% 	type	: type of B-splines
%		= '+ortho' (causal orthonormal, default), '-ortho' (anticausal  
%		orthonormal) or '*ortho' (symmetric orthonormal);
% 		= '+bspline' (causal B-spline), '-bspline' (anticausal 
%		B-spline) or '*bspline' (symmetric B-spline). 
% 		= '+dual' (causal dual), '-dual' (anticausal dual) or '*dual' 
% 		(symmetric dual). The last option is the flipped version of the
%		B-spline one.
% 	Output:
% 	FFTanalysisfilters	= [lowpassfilter;highpassfilter]	: FFT filter arrays
% 	FFTsynthesisfilters	= [lowpassfilter;highpassfilter] 	: FFT filter arrays
% 
% 	See also FFTWAVELETANALYSIS, FFTWAVELETSYNTHESIS.
% 	Uses FRACTSPLINEAUTOCORR
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

if alpha<=-0.5
	disp('The autocorrelation of the fractional splines exists only ')
	disp('for degrees strictly larger than -0.5!')
	FFTanalysisfilters=[];
	FFTsynthesisfilters=[];
	return
end
if M~=2^round(log(M)/log(2))
	disp(' ')
	disp('The size of the FFT must be a power of two!')
	disp(' ')
	FFTanalysisfilters=[];
	FFTsynthesisfilters=[];
	return
end
nu=0:1/M:(1-1/M);

if nargin<3
	type='+ortho';
end
if length(type)<=1
	error(['''' type '''' ' is an unknown filter type!'])
end

A=fractsplineautocorr(alpha,nu);
A2=[A A];
A2=A2(1:2:length(A2));		% A2(z) = A(z^2)

if type(2)=='o'|type(2)=='O'
	% orthonormal spline filters
	if type(1)=='*'
		lowa=sqrt(2)*abs((1+exp(-2*i*pi*nu))/2).^(alpha+1).*sqrt(A./A2);
	else
		if type(1)=='+'
			lowa=sqrt(2)*((1+exp(-2*i*pi*nu))/2).^(alpha+1).*sqrt(A./A2);
		else
			if type(1)=='-'
				lowa=sqrt(2)*((1+exp(2*i*pi*nu))/2).^(alpha+1).*sqrt(A./A2);
			else
				error(['''' type(1) '''' ' is an unknown filter prefix!'])
			end
		end
	end
	higha=exp(-2*i*pi*nu).*lowa;
	higha=conj([higha(M/2+(1:M/2)) higha(1:M/2)]);
	
	lows=lowa;
	highs=higha;
	FFTanalysisfilters=[lowa;higha];
	FFTsynthesisfilters=[lows;highs];
else
	% semi-orthonormal spline filters
	if type(1)=='*'
		lowa=sqrt(2)*abs((1+exp(-2*i*pi*nu))/2).^(alpha+1);
	else
		if type(1)=='+'
			lowa=sqrt(2)*((1+exp(-2*i*pi*nu))/2).^(alpha+1);
		else
			if type(1)=='-'
				lowa=sqrt(2)*((1+exp(2*i*pi*nu))/2).^(alpha+1);
			else
				error(['''' type '''' ' is an unknown filter type!'])
			end
		end
	end
	higha=exp(-2*i*pi*nu).*lowa.*A;
	higha=conj([higha(M/2+(1:M/2)) higha(1:M/2)]);
	
	lows=lowa.*A./A2;
	highs=higha./(A2.*[A(M/2+(1:M/2)) A(1:M/2)]);
	if type(2)=='d'|type(2)=='D' 
		% dual spline wavelets
		FFTanalysisfilters=[lowa;higha];
		FFTsynthesisfilters=[lows;highs];
	else
		% B-spline wavelets
		if type(2)=='b'|type(2)=='B'
			FFTsynthesisfilters=[lowa;higha];
			FFTanalysisfilters=[lows;highs];
		else
			error(['''' type '''' ' is an unknown filter type!'])
		end
	end	
end


