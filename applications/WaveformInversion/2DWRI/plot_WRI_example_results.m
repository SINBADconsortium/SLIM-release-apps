close all

%% set folders which contain results and select which results to load
startup
addpath([slimapps '/tools/functions/misc']);

%folder which contains results
current_dir=pwd;
resultsdir = '/results/';

% .mat files which contain results
result = 'workspace_workspaceWRI_BG_5_9HZ_2sweep';

load([current_dir resultsdir result]);

%% plot models

n=model.n;
x=model.x;
z=model.z;
v_est = reshape(1e3./sqrt(mk),model.n);
v0 = reshape(v0,n);

figure(1);
set(gca,'Fontsize',10)
imagesc(x,z,v_true,[1500 4500]);title('True velocity model')
xlabel('x [m]');ylabel('z [m]');axis equal tight;
colorbar;colormap jet

figure(2);
set(gca,'Fontsize',10)
imagesc(x,z,v0,[1500 4500]);title(['Initial velocity model'])
xlabel('x [m]');ylabel('z [m]');axis equal  tight
colorbar;colormap jet

figure(3);
set(gca,'Fontsize',10)
imagesc(x,z,(v_est),[1500 4500]);title(['Result WRI, \lambda=',num2str(params.lambda)])
xlabel('x [m]');colorbar;colormap jet
axis equal tight;


%% plot total objective function value (not normalized)
figure;
semilogy(cell2mat(f.obj'));pbaspect([4 1 1]);xlabel('iteration number');ylabel('objective value')
title('total objective function value');

%% First update at first frequency batch
v1 = reshape(1e3./sqrt(models{1}(:,2)),models_struct{1}.n);
v2 = reshape(1e3./sqrt(models{1}(:,3)),models_struct{1}.n);

figure(7);
set(gca,'Fontsize',10)
imagesc(x,z,v2-v1,[-200 200]);title('First update')
xlabel('x [m]');ylabel('z [m]');axis equal tight;

%% movie of the model estimate at every iteration
%some lines are commented out, which are for saving the movie to a file.
%These lines are partially dependent on the available movie formats (varies per matlab version and operating system) and may
%need to be adjusted.

models=cell2mat(models);

%writerObj = VideoWriter('Mov_wri_result','MPEG-4'); 
%writerObj.FrameRate = 10;
%open(writerObj);
counter=1;
figure
for i=1:size(models,2)
    
    err_rel(i)=norm((1e3./sqrt(models(:,i)))-v_true(:))/norm(v_true(:));
    if i>1
        evol(i)=norm(models(:,i)-models(:,i-1));
    end
    imagesc(model.x,model.z,reshape(1e3./sqrt(models(:,i)),model.n),[1500 5750]);colormap jet;;axis image
     pause(0.1)
    %Mov_frame(i) = getframe(gcf);
    %writeVideo(writerObj,Mov_frame(i));
    error_rel_wri(counter)=norm(1e3./sqrt(models(:,i))-1e3./sqrt(mk(:)))/norm(1e3./sqrt(mk(:)));
    counter=counter+1;
end
        
%% (normalized) difference with the true model
figure;set(gca,'Fontsize',10)
plot(nonzeros(error_rel_wri),'r','LineWidth',3);ylabel('Normalized error');
title('|| m estimate - m true ||_2');xlabel('Iteration nr.');pbaspect([3 1 1])

        
%% Crossections
%to be able to plot the results(estimated, startmodel and true model in one
%image, we need not resmaple them to a common grid):
[Xq Yq]=meshgrid(x,z);

[X Y]=meshgrid(model.x,model.z);
v_est_resampled = interp2(X,Y,v_est,Xq,Yq);

[X Y]=meshgrid(model.x,model.z);
v0_resampled = interp2(X,Y,v0,Xq,Yq);

l1=round(size(v_true,2)/4);
l2=round(size(v_true,2)/2); 
l3=round(2*size(v_true,2)/3);
figure;set(gca,'Fontsize',10)

subplot(1,3,1);set(gca,'Fontsize',10);plot(v_true(:,l1),z,'LineWidth',1.5);set(gca,'YDir','reverse');hold on;title(['x =',num2str(x(l1)),'[m]']);xlabel('Velocity [m/s]');ylabel('z [m]');axis([1000 5000 0  max(z)+100]);pbaspect([0.5 1 1])
subplot(1,3,2);set(gca,'Fontsize',10);plot(v_true(:,l2),z,'LineWidth',1.5);set(gca,'YDir','reverse');hold on;title(['x =',num2str(x(l2)),'[m]']);xlabel('Velocity [m/s]');ylabel('z [m]');axis([1000 5000 0  max(z)+100]);pbaspect([0.5 1 1])
subplot(1,3,3);set(gca,'Fontsize',10);plot(v_true(:,l3),z,'LineWidth',1.5);set(gca,'YDir','reverse');hold on;title(['x =',num2str(x(l3)),'[m]']);xlabel('Velocity [m/s]');ylabel('z [m]');axis([1000 5000 0  max(z)+100]);pbaspect([0.5 1 1])

subplot(1,3,1);plot(v_est_resampled(:,l1),z,'r','LineWidth',1.5);set(gca,'YDir','reverse');hold on;
subplot(1,3,2);plot(v_est_resampled(:,l2),z,'r','LineWidth',1.5);set(gca,'YDir','reverse');hold on;
subplot(1,3,3);plot(v_est_resampled(:,l3),z,'r','LineWidth',1.5);set(gca,'YDir','reverse');hold on;

hold on;
subplot(1,3,1);plot(v0_resampled(:,l1),z,'g','LineWidth',1.5);set(gca,'YDir','reverse');hold on;
subplot(1,3,2);plot(v0_resampled(:,l2),z,'g','LineWidth',1.5);set(gca,'YDir','reverse');hold on;
subplot(1,3,3);plot(v0_resampled(:,l3),z,'g','LineWidth',1.5);set(gca,'YDir','reverse');hold on;

p=legend('True','estimate','start','Location','Southoutside');

