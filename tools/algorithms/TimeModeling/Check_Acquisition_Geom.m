function [] = Check_Acquisition_Geom(model1)
    %% Checks acquisition geometry :
    %  Receiver in the domain
    %  Source in domain
    %  Source and receiver arrays have the correct size
	% Author : Mathias Louboutin
	%         Seismic Laboratory for Imaging and Modeling
	%         Department of Earth, Ocean, and Atmosperic Sciences
	%         The University of British Columbia
    % Date : July, 2015
    % Check receiver position
    minrecx=min(model1.xrec(:));
    maxrecx=max(model1.xrec(:));

    if model1.n(3)>1
        minrecy=min(model1.yrec(:));
        maxrecy=max(model1.yrec(:));
        range= (minrecx>=model1.o(2)) && (minrecy>=model1.o(3)) && ...
               (maxrecx<=(model1.o(2)+(model1.n(2)-1)*model1.d(2))) && ...
               (maxrecy<=(model1.o(3)+(model1.n(3)-1)*model1.d(3))); 
    else
         range= (minrecx>=model1.o(2))&& ...
               (maxrecx<=(model1.o(2)+(model1.n(2)-1)*model1.d(2)));
    end

    assert(range,'Receiver outside the computational domain, make sure the model fits the acquisition setup');

    % Check source position

    minsrcx=min(model1.xsrc(:));
    maxsrcx=max(model1.xsrc(:));

    if model1.n(3)>1
        minsrcy=min(model1.ysrc(:));
        maxsrcy=max(model1.ysrc(:));
        range= (minsrcx>=model1.o(2)) && (minsrcy>=model1.o(3)) && ...
               (maxsrcx<=(model1.o(2)+(model1.n(2)-1)*model1.d(2))) && ...
               (maxsrcy<=(model1.o(3)+(model1.n(3)-1)*model1.d(3))); 
    else
         range= (minsrcx>=model1.o(2))&& ...
               (maxsrcx<=(model1.o(2)+(model1.n(2)-1)*model1.d(2)));
    end
    
    assert(range,'Source outside the computational domain, make sure the model fits the acquisition setup');

    % Check receiver size

    if strcmp(model1.type,'marine')
        assert(size(model1.xrec,2)==size(model1.xsrc,2) || size(model1.yrec,2)==size(model1.ysrc,2),...
        ['For marine acquisition, receiver position array must be nsrc*nrec' ...
        'You can fill the receiver array with the last value to force every shot to have the same number of receivers']);
    end

    assert(size(model1.xsrc,1)==size(model1.ysrc,1) || size(model1.xsrc,2)==size(model1.ysrc,2) || ...
    size(model1.xsrc,1)==size(model1.zsrc,1) || size(model1.xsrc,2)==size(model1.zsrc,2),'X,Y and Z source arrays must have the same size');


    if isfield(model1,'tfire')
        assert(length(model1.tfire)==length(model1.xsrc(:)),'Number of firing times must be equal to the number of source positions for continuous acquisition');
    end

end

