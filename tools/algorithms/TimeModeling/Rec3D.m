function P = Rec3D(model)
	%% 3D projection of the wavefield onto the receiver location
    %  Receiver data is read from a 3x3x3 cube with exponential damping
	% Author : Mathias Louboutin
	%         Seismic Laboratory for Imaging and Modeling
	%         Department of Earth, Ocean, and Atmosperic Sciences
	%         The University of British Columbia
    % Date : July, 2015
	nrec=length(model.rlx_loc(:));

	if nrec==0
		P=[];
	else
	
	
		[uniz,iz]=unique(model.rlz_loc(:));
		[unixx,ix]=unique(model.rlx_loc(:));
		[uniy,iy]=unique(model.rly_loc(:));

		if model.n(3)>1

            Pz=exp(-1)*getLA(model.mmz_loc,uniz-model.d(1))+getLA(model.mmz_loc,uniz)+exp(-1)*getLA(model.mmz_loc,uniz+model.d(1));
            Px=exp(-1)*getLA(model.mmx_loc,unixx-model.d(2))+getLA(model.mmx_loc,unixx)+exp(-1)*getLA(model.mmx_loc,unixx+model.d(2));
            Py=exp(-1)*getLA(model.mmy_loc,uniy-model.d(3))+getLA(model.mmy_loc,uniy)+exp(-1)*getLA(model.mmy_loc,uniy+model.d(3));
            
            P=kron(Py,kron(Px,Pz));
            P=sparse(P);
		else
            Pz=exp(-1)*getLA(model.mmz_loc,uniz-model.d(1))+getLA(model.mmz_loc,uniz)+exp(-1)*getLA(model.mmz_loc,uniz+model.d(1));                      
            Px=exp(-1)*getLA(model.mmx_loc,unixx-model.d(2))+getLA(model.mmx_loc,unixx)+exp(-1)*getLA(model.mmx_loc,unixx+model.d(2));

            P=kron(Px,Pz);
            P=sparse(P);
        end
	end
end