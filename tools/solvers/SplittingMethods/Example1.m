% visual example to see how alternating projection methods work (3 - sets)
% Bas Peters, Feb. 2016
close all
clear
 
%% define constraints

%constraint set 1: (x+y)^2< or equal to r
r=sqrt(9);
P1 = @(input) input/norm(input,2) * min(r,norm(input,2));

%constraint set 2: x< or equal to 2 and y<or equal to 3
d1=3;
d2=2;
P2 = @(input) vec(median([[d1;d2] input [-d1;-d2]]'));

%constraint set 3: projection onto a hyperplane
A(1,1)=1;
A(1,2)=-1;

xt=[1.9;2.1];
b=0;

P3 = @(input)  input - ((A(1,:)*input-b(1))/norm(A(1,:))^2).*A(1,:)';% + randn(1)*0.1;

Proj{1}=P1;
Proj{2}=P2;
Proj{3}=P3;


%% initial guess
v=[2.0 ; 2.5];

%% parallel proximal dykstra-like splitting
x_par_dyk(:,1)=v;
options_dyk.log_vec=1;
counter=2;
for i=1:20
options_dyk.maxIt=i;
 [x_par_dyk(:,counter),res,evol_rel,nx,np,nz,xlog,zlog,plog]=Dykstra_prox_parallel(v,Proj,options_dyk);
 counter=counter+1;
end
x_par_dyk=x_par_dyk';
err_d_par=x_par_dyk-repmat([2 2],size(x_par_dyk,1),1);


%% Plot solution path
sq=[-3  2;
     3  2;
     3 -2;
    -3 -2;];

figure(1);set(gca,'Fontsize',14)
plot(sq(:,1),sq(:,2),'b','LineWidth',3);hold on

h=ezplot('x.^2 + y.^2=9',[1 3]);set(h, 'Color', 'b','LineWidth',2);hold on
h=ezplot('1*x-1*y=0',[1 3]);set(h, 'Color', 'k','LineWidth',2);hold on
title('')

plot(x_par_dyk(:,1),x_par_dyk(:,2),'-.xm','LineWidth',1.5);hold on
axis equal
axis([1.8 2.5 1.8 2.6])

%% plot residual and error
figure;
semilogy(sqrt(sum(err_d_par.^2,2)));hold on
semilogy(res,'r');
xlabel('iteration number');title('Parallel proximal Dykstra')
legend('error','residual')
