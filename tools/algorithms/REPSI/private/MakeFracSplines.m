function [w_analysis,w_synthesis]= MakeFracSplines(J,bandnumber,alpha,type);
%
%INPUT:
% J             - number of octaves
% bandnumber    - octave #
% alpha         - vector with orders alpha>-0.5 : degree of the fractional splines 
% 	type	- type of B-splines
%		= '+ortho' (causal orthonormal, default), '-ortho' (anticausal  
%		orthonormal) or '*ortho' (symmetric orthonormal);
% 		= '+bspline' (causal B-spline), '-bspline' (anticausal 
%		B-spline) or '*bspline' (symmetric B-spline). 
% 		= '+dual' (causal dual), '-dual' (anticausal dual) or '*dual' 
% 		(symmetric dual). The last option is the flipped version of the
%		B-spline one.
%
%OUTPUT:
% w             - matrix with the wavelets
if min(alpha)<=-0.5
	disp('The autocorrelation of the fractional splines exists only ')
	disp('for degrees strictly larger than -0.5!')
	w=[];
	return
end
if nargin<4
	type   =   '+ortho';
end

nalpha         =   length(alpha);
y              =   zeros(1,2^J);
w_analysis     =   zeros(2^J,nalpha);
w_synthesis    =   zeros(2^J,nalpha);
for ialpha=1:nalpha
  buf           =  dirac(2^(J-bandnumber),1);
  wbuf          =  winsert(y,buf',J,bandnumber);
  [FFTanalysisfilters,FFTsynthesisfilters]=FFTfractsplinefilters(2^J,alpha(ialpha), type);
  wbuf1         =  FFTwaveletsynthesis(wbuf,FFTsynthesisfilters,J);
  wbuf2         =  FFTwaveletsynthesis(wbuf,FFTsynthesisfilters,J);
  w_analysis(:,ialpha) = wbuf1';
  w_synthesis(:,ialpha) = wbuf2';
end
