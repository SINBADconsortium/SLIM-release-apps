function [qmf] = MakeFracSplinesFilters(M,alpha,type)

%-----------------------------------------------------------------------------------
% INPUT:
% M             - length of the filters
% alpha         - vector with orders alpha>-0.5 : degree of the fractional splines 
% type	        - type of B-splines
%		  '+ortho'   (causal orthonormal, default)
%                 '-ortho'   (anticausal  orthonormal)  
%                 '*ortho'   (symmetric orthonormal);
% 		  '+bspline' (causal B-spline)
%                 '-bspline' (anticausal B-spline)
%                 '*bspline' (symmetric B-spline). 
% 		  '+dual'    (causal dual)
%                 '-dual'    (anticausal dual)
%                 '*dual'    (symmetric dual). 
%
% The last option is the flipped version of the	B-spline one.
%
% OUTPUT:
% qmf           - matrix with the wavelets
%
% USAGE:
%          [qmf] = MakeFracSplinesFilters(M,alpha,type);  
%
%   	References:
% 	[1] M. Unser and T. Blu, "Fractional splines and wavelets," 
% 	SIAM Review, in press.
% 	[2] M. Unser and T. Blu, "Construction of fractional spline wavelet bases," 
% 	Proc. SPIE vol 3813, Wavelet Applications in Signal and Image 
% 	Processing VII, in press. 
%-----------------------------------------------------------------------------------

if min(alpha)<=-0.5
	disp('The autocorrelation of the fractional splines exists only ')
	disp('for degrees strictly larger than -0.5!')
	w=[];
	return
end
if nargin<3
	type   =   '+ortho';
end

% check type
if iscell(type),
  I            =   length(type);
  typ          =   type;
else
  I            =   1;
  typ          =   cell(I);
  typ{1}       =   type;
end
N              =   length(alpha);
qmf            =   zeros(N*I,M);
%qmf            =   zeros(N*I*2,M);

for n=1:N
  for i=1:I
    [F        ]=FFTfractsplinefilters(M,alpha(n), ...
						  typ{i});
    idx        =  (n-1)*I+(i-1)+1;
    qmf(idx,:)  =  real(ifft(F(1,:)));    
    if ((alpha(n)>=0 & alpha(n)<2)),
      switch typ{i}
       case '*ortho'
	shif = 11;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '-ortho'
	shif = 12;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '+ortho'
	shif = 10;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
      end
    end
    if ((alpha(n)>=2 & alpha(n)<3)),
      switch typ{i}
       case '*ortho'
	shif = 11;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '-ortho'
	shif = 13;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '+ortho'
	shif = 9;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
      end
    end
    if ((alpha(n)>=3 & alpha(n)<4)),
      switch typ{i}
       case '*ortho'
	shif = 20;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '-ortho'
	shif = 22;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '+ortho'
	shif = 18;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
      end
    end
    if ((alpha(n)>=4 & alpha(n)<5)),
      switch typ{i}
       case '*ortho'
	shif = 20;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '-ortho'
	shif = 23;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '+ortho'
	shif = 17;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
      end
    end
    if ((alpha(n)>=5 & alpha(n)<6)),
      switch typ{i}
       case '*ortho'
	shif = 29;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '-ortho'
	shif = 32;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '+ortho'
	shif = 26;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
      end
    end
    if ((alpha(n)>=6 & alpha(n)<7)),
      switch typ{i}
       case '*ortho'
	shif = 29;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '-ortho'
	shif = 33;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
       case '+ortho'
	shif = 25;
	qmf(idx,:)  = shift(qmf(idx,:)',shif)';
      end
    end
  end
end
