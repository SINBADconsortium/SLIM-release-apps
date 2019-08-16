function [alpha] = WolfeLS(fh,x0,f0,g0,alpha0,Z)
%% Line search for weak Wolf conditions
% (http://cs.nyu.edu/overton/mstheses/skajaa/msthesis.pdf, algorihtm 3)
%
% Author : Tristan van Leuwen edited by Mathias Louboutin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date : July, 2015
if nargin<6
	Z=1;
end
fprintf(1,'  # iter, stepsize   , ft         , gt*s0    ,f0+c1*alpha*g0*s0,  c2 g0*s0\n');
%fprintf(1,'  # iter, stepsize   , f(x)       , f0       , ||g0||_2\n');
lsiter = 0;
c1 = 1e-2;
c2 = 0.9;
done = 0;
mu = 0;
nu = inf;
alpha =alpha0;
alphap=0;
s0=-Z*g0;
while ~done
    if nu < inf
        alpha = (nu + mu)/2;
    else
        alpha = 2*alpha;
    end
    
    if lsiter < 10
        [ft,gt] = fh(x0+alpha*s0);
        lsiter = lsiter + 1;
    else
		alpha=alphap;
        break;
    end
    
	fprintf(1,'      >%d, %1.5e, %1.5e, %1.5e, %1.5e, %1.5e\n',lsiter, alpha, ft, gt(:)'*s0(:),f0 + c1*alpha*g0(:)'*s0(:),c2*g0(:)'*s0(:));
    %fprintf(1,'      >%d, %1.5e, %1.5e, %1.5e, %1.5e\n',lsiter, alpha, ft,f0,norm(g0));
    
    if ft > f0 + c1*alpha*g0(:)'*s0(:)
        nu = alpha;
    elseif gt(:)'*s0(:) < c2*g0(:)'*s0(:)
        mu = alpha;
		alphap=alpha;
    else
		if(isnan(ft) || isinf(ft))
			if lsiter==10
				alpha=alphap;
				done=1;
			else
				nu=alpha;
			end
		else
			done=1;
		end
    end
	
end

end
