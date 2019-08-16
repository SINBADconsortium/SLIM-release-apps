%This is the main script for using scaled gradient projection (SGP) to 
%minimize the wavefield reconstruction inversion (WRI) objective
%subject to total variation and spatially varying bound constraints.
%It minimizes data misfit plus pde misfit plus subject to model constraints
%min_{m,u} .5*||Pu-d||^2 + .5*lam^2||A(m)u-q||^2 
%s.t. ||m||_TV <= tau and m_ij in [b_ij,B_ij].
%This extends Tristan van Leeuwen's earlier penalty method implementation.

%setup model and parameters for a particular example
example = 'blob'; %simple model with high velocity blob and smooth background
%example = 'BPC'; %example using downsampled center part of BP 2004 Salt Structure Data
[model,pm,Ps,Pr,q,d,ssW] = setup_problem(example);

%preliminary plots
figure(51); clf; 
imagesc(model.xt,model.zt,reshape(sum(Ps,1)' + 2*sum(Pr,1)',model.n)); 
title('source and receiver locations'); figadjust

%loop over frequency batches, each time calling pTV to apply SGP
%store results in mb
mb = zeros(model.n(1),model.n(2),model.batches);
minit = model.minit; 
%project onto feasible set if outside
if (pm.TV(minit) > pm.tau || max((minit(:)-model.mmin(:)).*(minit(:)-model.mmax(:))) > 0)
    disp('projecting initial model onto constraint set');
    minit = proj_TVbounds(minit,.99*pm.tau,model.h,pm.dpw,model.mmin,model.mmax);
end
mb(:,:,1) = minit; 
cinit = 1e-12; %initial damping parameter (will be automatically updated)
for b = 1:model.batches
    %solve problem for current frequency batch, using previous as init
    [mTV,energy,oits,cinit] = pTV(model,pm,b,mb(:,:,b),Ps,Pr,q,d,ssW,cinit);
    mb(:,:,b+1) = mTV; %store intermediate results using different subsets of frequencies
    
    %plot objective for latest frequency batch
    figure(71)
    clf;
    plot(energy); title('objective versus iteration');
    
    %show intermediate velocity estimate
    figure(72)
    clf;
    imagesc(model.xt,model.zt,1./sqrt(mTV)); 
    title(['velocity estimate, ' num2str(oits) ' outer iterations, ' num2str(b) ' freq batches']);
    colorbar; caxis([min(model.vmin(:)) max(model.vmax(:))]);
    drawnow;
end

%show final velocity estimate
figure(73); clf;
imagesc(model.xt,model.zt,1./sqrt(mTV)); title('final velocity estimate');
colorbar; caxis([min(model.vmin(:)) max(model.vmax(:))]);

%compare to true velocity
figure(74); clf;
imagesc(model.xt,model.zt,1./sqrt(model.mtrue)); title('true velocity');
colorbar; caxis([min(model.vmin(:)) max(model.vmax(:))]);

%and initial velocity
figure(75); clf;
imagesc(model.xt,model.zt,1./sqrt(model.minit)); title('initial velocity');
colorbar; caxis([min(model.vmin(:)) max(model.vmax(:))]);

%plot relative error of model (normalized to start at one)
figure(76); clf;
model_error = sqrt(sum((reshape(mb,model.N,model.batches+1) - ...
    repmat(model.mtrue(:),1,model.batches+1)).^2,1))/...
    norm(minit(:) - model.mtrue(:));
plot(0:model.batches,model_error); title('relative model error vs frequency batch');

%make a movie of the results with one frame for each frequency batch
%make_movie %script to make movie from mb
%movie(99,velocity_movie,1,4); %play movie at 4 fps

%save results...
    
