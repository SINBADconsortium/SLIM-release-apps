function [l p] = plot_MeanCmp(x,m1,std1,m2,std2)

%[l,p]=boundedline(m1,x,std1,'-*', m2, x, std2, '-or', 'orientation','horiz','alpha');
B(:,:,1) = [std1(:),std1(:)];
B(:,:,2) = [std2(:),std2(:)];
[l,p]=boundedline([m1(:),m2(:)],x,  B,'alpha','-','orientation','horiz');
set(l,'linewidth',2)
set(p,'linewidth',2)
outlinebounds(l,p);
ylim([x(1) x(end)]);
xlim([0,6]);
xlabel('Velocity [km/s]');
ylabel('Depth [m]');
set(gca,'Ydir','reverse')
set(gca,'fontsize',20)

