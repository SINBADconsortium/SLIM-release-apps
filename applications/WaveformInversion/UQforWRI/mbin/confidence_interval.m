function y = confidence_interval(x,alpha);
	%% y = confidence_interval(x,alpha);
	% Generate the alpha confidence interval for the vector x
	% 
	% Input: 
	% x     - sample vector
	% alpha - Confidence leval 
	% 
	% Output
	% y     - confidence interval s.t. P(x < y(1)) = (1 - alpha)/2; P(x>y(2)) = (1 - alpha)/2;
	% 
	% Written by Zhilong Fang, SLIM, UBC, 2013/11/15
	
	beta = (1-alpha)/2;
	nout = floor(length(x)*beta)+1;
	z    = sort(x);
	y    = [z(nout+1),z(length(x)-nout)];
	 
	 
	 
	 
	 
	 