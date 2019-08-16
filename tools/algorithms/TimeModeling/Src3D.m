function P = Src3D(model,multi)
	%% 3D projection of the source onto the computational grid
    %  Source is injected on a 3x3x3 cube with exponential damping to avoid
    %  dispersion from a point source (dirac)
	% Author : Mathias Louboutin
	%         Seismic Laboratory for Imaging and Modeling
	%         Department of Earth, Ocean, and Atmosperic Sciences
	%         The University of British Columbia
    % Date : July, 2015
	N=length(model.mmx_loc)*length(model.mmy_loc)*length(model.mmz_loc);
	nrec=length(model.slx_loc(:));
	if nargin==1
		multi=0;
	end
	if nrec==0
		P=[];
    else
    d=1.5*model.dt;

	[uniz,iz]=unique(model.slz_loc(:));
	[unixx,ix]=unique(model.slx_loc(:));
	[uniy,iy]=unique(model.sly_loc(:));


	P=sparse(1,N);
	if model.n(3)>1
        Pz=1/d*getLA(model.mmz_loc,uniz-d)+getLA(model.mmz_loc,uniz)+1/d*getLA(model.mmz_loc,uniz+d);
        Px=1/d*getLA(model.mmx_loc,unixx-d)+getLA(model.mmx_loc,unixx)+1/d*getLA(model.mmx_loc,unixx+d);
        Py=1/d*getLA(model.mmy_loc,uniy-d)+getLA(model.mmy_loc,uniy)+1/d*getLA(model.mmy_loc,uniy+d);

        P=kron(Py,kron(Px,Pz));
        P=sparse(P);
	else
        Pz=1/d*getLA(model.mmz_loc,uniz-d)+getLA(model.mmz_loc,uniz)+1/d*getLA(model.mmz_loc,uniz+d);
        Px=1/d*getLA(model.mmx_loc,unixx-d)+getLA(model.mmx_loc,unixx)+1/d*getLA(model.mmx_loc,unixx+d); 

        P=kron(Px,Pz);
        P=sparse(P);
	end
		%Computational source grid interpolated from the physical source grid
	end
end
