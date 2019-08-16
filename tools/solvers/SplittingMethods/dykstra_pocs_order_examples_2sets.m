%Generate visual example about different projection methods (Dykstra vs
%POCS). Also shows what the effect of set ordering is.
% Bas Peters, 2016


close all
clear

%constraint set 1: (x+y)^2< or equal to r
%Projector 1:
r=sqrt(9);
P1 = @(input) input/norm(input,2) * min(r,norm(input,2));

%constraint set 2: x< or equal to 2 and y<or equal to 3
%Projector 1:
d1=3;
d2=2;
P2 = @(input) vec(median([[d1;d2] input [-d1;-d2]]'));


%initial guess
v=[2.5 ; 3.0];
%v=[2.0 ; 2.5];


%% POCS (ordering 1)
x_pocs(:,1)=v;
counter=2;
for i=1:10
    x_pocs(:,counter)=P1(x_pocs(:,counter-1));
    counter=counter+1;
    x_pocs(:,counter)=P2(x_pocs(:,counter-1));
    counter=counter+1;
end

%% POCS (ordering 2)

x_pocs2(:,1)=v;
counter=2;
for i=1:10
    x_pocs2(:,counter)=P2(x_pocs2(:,counter-1));
    counter=counter+1;
    x_pocs2(:,counter)=P1(x_pocs2(:,counter-1));
    counter=counter+1;
end

%% dykstra (ordering 1)
out(:,1)=v;
counter=2;
x=v;
p=[0;0];
q=[0;0];
for i=2:10
options_dyk.maxIt=i;
y=P1(x+p);
out(:,counter)=y;
counter=counter+1;
p=x+p-y;
x=P2(y+q);
out(:,counter)=x;
counter=counter+1;
q=y+q-x;
end
out=out';

%% dykstra (ordering 2)

out2(:,1)=v;
counter=2;
x=v;
p=[0;0];
q=[0;0];
for i=2:10
options_dyk.maxIt=i;
y=P2(x+p);
out2(:,counter)=y;
counter=counter+1;
p=x+p-y;
x=P1(y+q);
out2(:,counter)=x;
counter=counter+1;
q=y+q-x;
end
out2=out2';

%% plot results
sq=[-3  2;
     3  2;
     3 -2;
    -3 -2;];
    
%plot this system and its solution
figure(1);
subplot(1,2,1); set(gca,'Fontsize',14)
plot(x_pocs(1,:),x_pocs(2,:),'-ok','LineWidth',3.5);hold on;
plot(out(:,1),out(:,2),'-.^r','LineWidth',1.5)
plot(sq(:,1),sq(:,2),'b','LineWidth',3);hold on
h=ezplot('x.^2 + y.^2=9',[1 3]);set(h, 'Color', 'b','LineWidth',2);hold on
title('')
legend('POCS 1','Dykstra 1','Location','NorthWest')
axis([1.8 2.6 1.8 3.1])
%labels = cellstr( num2str([1:20]') );  %' # labels correspond to their order
%text(out(:,1),out(:,2),labels,'VerticalAlignment','bottom', ...
%                             'HorizontalAlignment','right', 'FontSize',20)
subplot(1,2,2); set(gca,'Fontsize',14)
plot(x_pocs2(1,:),x_pocs2(2,:),'-ok','LineWidth',3.5);hold on;
plot(out2(:,1),out2(:,2),'-.^r','LineWidth',1.5)
plot(sq(:,1),sq(:,2),'b','LineWidth',3);hold on
h=ezplot('x.^2 + y.^2=9',[1 3]);set(h, 'Color', 'b','LineWidth',2);hold on
title('')
legend('POCS 2','Dykstra 2','Location','NorthWest');
axis([1.8 2.6 1.8 3.1])




%export_fig('Dykstra_vs_POCS_ordering','-pdf','-transparent','-depsc');%close all
