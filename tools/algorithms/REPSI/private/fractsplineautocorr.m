function A=fractsplineautocorr(alpha,nu) 

% FRACTSPLINEAUTOCORR Frequency domain computation of fractional spline 
% 	autocorrelation.  A=fractsplineautocorr(alpha,nu) computes the 
% 	frequency response of the autocorrelation filter A(exp(2*i*Pi*nu)) 
% 	of a fractional spline of degree alpha.  It uses an acceleration 
% 	technique to improve the convergence of the infinite sum by 4 
% 	orders.
% 
% 	See also FFTSPLINEFILTERS
% 	
% 	Author: Thierry Blu, October 1999 
% 	Biomedical Imaging Group, EPFL, Lausanne, Switzerland.  
% 	This software is downloadable at http://bigwww.epfl.ch/
% 	
% 	References:
% 	[1] M. Unser and T. Blu, "Fractional splines and wavelets," 
% 	SIAM Review, in press.
	
N=100;			% number of terms of the summation for computing
				% the autocorrelation frequency response

if alpha<=-0.5
	disp('The autocorrelation of the fractional splines exists only ')
	disp('for degrees strictly larger than -0.5!')
	A=[];
	return
end
				
S=zeros(1,length(nu));
err=[];
err0=[];
for n=-N:N
	S=S+abs(sinc(nu+n)).^(2*alpha+2);
end
U=2/(2*alpha+1)/N^(2*alpha+1);
U=U-1/N^(2*alpha+2);
U=U+(alpha+1)*(1/3+2*nu.*nu)/N^(2*alpha+3);
U=U-(alpha+1)*(2*alpha+3)*nu.*nu/N^(2*alpha+4);
U=U.*abs(sin(pi*nu)/pi).^(2*alpha+2);
A=S+U;
