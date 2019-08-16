vt     = [2000; 2500];
dv     = [10];
nm1    = 10;
nm2    = 8;
PaperSize = [8,6];
Fsize     = 16;

AZ = -77;
EL = 10;




for i = 1:length(dv)
    dvt   = dv(i);
    v_1 = [vt(1)-dvt*nm1:dvt:vt(1)+dvt*nm1];
    v_2 = [vt(2)-dvt*nm1:dvt:vt(2)+dvt*nm2];
    [x_2 x_1] = meshgrid(v_2,v_1);
    f_str = ['WRIvsFWI_' num2str(dvt) '.mat'];
    load(f_str);
    figure;h(1)=surf(x_2,x_1,f_fwi);set(h(1),'FaceColor','red');
   % hold on; h(2) = surf(x_2,x_1,f_fwi_app);set(h(2),'FaceColor','blue');
    %legend('True misfit', 'Quadratic approximation','Location','North');
    view(AZ,EL);set(gca,'fontsize', Fsize);
    xlabel('\bf{v}_{2} [m/s]','fontsize', Fsize);
    ylabel('\bf{v}_{1} [m/s]','fontsize', Fsize);
    zlabel('\bf{Normalized misfit value}','fontsize', Fsize);
    zlimt = zlim;
    xlim  = [v_2(1),v_2(end)];
    ylim  = [v_1(1),v_1(end)];
    fig_str = ['../fig/WRIvsFWI4/FWI_APP_dv_' num2str(dvt)];
    fid     = gcf;
    PrintFigure(fid,PaperSize,fig_str);
    pause(0.5)
    figure;h(1)=surf(x_2,x_1,f_wri);set(h(1),'FaceColor','red');
    %hold on; h(2) = surf(x_2,x_1,f_wri_app);set(h(2),'FaceColor','blue');
    %legend('True misfit', 'Quadratic approximation','Location','North');
    view(AZ,EL);set(gca,'fontsize', Fsize);
    xlabel('\bf{v}_{2} [m/s]','fontsize', Fsize);
    ylabel('\bf{v}_{1} [m/s]','fontsize', Fsize);
    zlabel('\bf{Normalized misfit value}','fontsize', Fsize);
    fig_str = ['../fig/WRIvsFWI4/WRI_APP_dv_' num2str(dvt)];
    zlim(zlimt);
    xlim  = [v_2(1),v_2(end)];
    ylim  = [v_1(1),v_1(end)];
    fid     = gcf;
    PrintFigure(fid,PaperSize,fig_str);
end

