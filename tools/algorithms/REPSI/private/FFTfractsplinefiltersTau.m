function [FFTanalysisfilters,FFTsynthesisfilters]= ...
    FFTfractsplinefiltersTau(M,alpha,tau,type)
%	Usage:
%
%	[FFTanalysisfilters,FFTsynthesisfilters]=FFTfractsplinefiltersTau(M,alpha,tau,type) 
% 	Provides the frequency response of the low- and high-pass 
% 	filters that generate the orthonormal or semi-orthonormal
% 	(B-spline  orr dual) fractional splines of degree alpha, shift
% 	tau and of given type.
% 	See also FFTWAVELETANALYSIS, FFTWAVELETSYNTHESIS.
% 	Uses FRACTSPLINEAUTOCORR%% 	Author: Thierry Blu, October
% 	1999, Revised January 2001
% 	Biomedical Imaging Group, EPFL, Lausanne, Switzerland.
% 	This software is downloadable at http://bigwww.epfl.ch/
%
% 	M 	: size of the input signal = length of FFTfilters = 2^N
% 	alpha 	: degree of the fractional splines, must be >-0.5 
% 	tau	: shift or asymmetry of the fractional splines, we 
%       suggest to restrict this value to the interval [-1/2,+1/2] 
%       because tau+1 leads to the same wavelet space as tau.
%	Particular cases
%		 	tau=0 <=> symmetric splines; 
%			tau=±1/2 <=> max. dissymetric splines
%   			tau=(alpha+1)/2 <=> causal splines)
% 	type	: type of the B-splines
%		= 'ortho' (orthonormal, default)
% 		= 'bspline' (B-spline) 
% 		= 'dual' (dual). The last option is the flipped version
% 		of the B-spline one.
%
%		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%		%%%%%%%%%%%% OUTPUT %%%%%%%%%%%	
%		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%
% 	FFTanalysisfilters	= [lowpassfilter;highpassfilter]	:
% 	FFT filter arrays
% 	FFTsynthesisfilters	= [lowpassfilter;highpassfilter] 	:
% 	FFT filter arrays
% 
%		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%		%%%%%%%%%%%% REFERENCES %%%%%%%
%		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%
% 	[1] M. Unser and T. Blu, "Fractional splines and wavelets," SIAM
% 	Review, vol. 42, no. 1, pp. 43-67, January 2000. 
% 	[2] M. Unser and T. Blu, "Construction of fractional spline
% 	wavelet bases," Proc. SPIE, Wavelet Applications in Signal and
% 	Image
%	Processing VII, Denver, CO, USA, 19-23 July, 1999,
%		vol. 3813, pp. 422-431.
% 	[3] T. Blu and M. Unser, "The fractional spline wavelet
% 	transform: definition and implementation," Proc. IEEE
% 	International
%		Conference on Acoustics, Speech, and Signal Processing
%		(ICASSP'2000), Istanbul, Turkey, 5-9 June 2000,
%		vol. I, pp. 512-515 . 

u=alpha/2-tau;v=alpha/2+tau;
if alpha<=-0.5	
  disp('The autocorrelation of the fractional splines exists only ')	
  disp('for degrees strictly larger than -0.5!')	
  FFTanalysisfilters=[];	
  FFTsynthesisfilters=[];	
  return
end
  if M~=2^round(log(M)/log(2))	
  disp(' ');
  disp('The size of the FFT must be a power of two!');	
  disp(' ');
  FFTanalysisfilters=[];	
  FFTsynthesisfilters=[];	
  return
end
nu=0:1/M:(1-1/M);
A=fractsplineautocorr(alpha,nu);
A2=[A A];A2=A2(1:2:length(A2));		
% A2(z) = A(z^2)
if type(1)=='o'|type(1)=='O'	% orthonormal spline filters	
  lowa=sqrt(2)*((1+exp(2*i*pi*nu))/2).^(u+0.5).*((1+exp(-2*i*pi*nu))/ ...
                                                 2).^(v+0.5).*sqrt(A./A2);
  higha=exp(2*i*pi*nu).*lowa;	higha=conj([higha(M/2+(1:M/2)) higha(1:M/ ...
                                                    2)]);
  lows=lowa;	
  highs=higha;	
  FFTanalysisfilters=[lowa;higha];	
  FFTsynthesisfilters=[lows;highs];
else	
% semi-orthonormal spline filters	

   lowa=sqrt(2)*((1+exp(2*i*pi*nu))/2).^(u+0.5).*((1+exp(-2*i*pi*nu))/2).^(v+0.5);
     
   higha=exp(2*i*pi*nu).*lowa.*A;	
   higha=conj([higha(M/2+(1:M/2)) higha(1:M/2)]);
   lows=lowa.*A./A2;	
   highs=higha./(A2.*[A(M/2+(1:M/2)) A(1:M/2)]);
   if type(1)=='d'|type(1)=='D' 		
% dual spline wavelets		
     FFTanalysisfilters=[lowa;higha];
     FFTsynthesisfilters=[lows;highs];	
   else		
% B-spline wavelets
      if type(1)=='b'|type(1)=='B'
        FFTsynthesisfilters=([lowa;higha]);
        FFTanalysisfilters=([lows;highs]);
      else			
        error(['''' type '''' ' is an unknown filter type!'])
      end
   end	
end