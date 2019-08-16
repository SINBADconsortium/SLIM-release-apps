function mycurtainfun(R,T)
%
m = size(R);
figure('WindowButtonDownFcn',@WinBtnDownCallback,'WindowButtonMotionFcn',@WinBtnMotionCallback0)


hax = gca;

[xd,yd]=ndgrid([0.5 m(1)+0.5],[0.5 m(2)+0.5]);
hR=surface(yd,xd,0*xd);
hT=surface(yd,xd,0*xd);
pos = round(m(2)/2);
C   = uint8(0*T); C(:,1:pos)=1;
set(hR,'facecolor','texture','facealpha','texturemap','cdata',double(R),'alphadata',C);
set(hT,'facecolor','texture','facealpha','texturemap','cdata',double(T),'alphadata',1-C);
hold on;
hcurtain(1) = plot([ m(2)/2 m(2)/2 ] , [ 0.5 m(1)+0.5 ] , 'r','linewidth',3);
hcurtain(2) = plot([ m(2)/2 m(2)/2 ] , [ 0.5 m(1)+0.5 ] , 'y','linewidth',1);
hold off;
axis ij tight;

myptr = [
    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
    NaN   NaN   NaN   NaN   2     2     NaN   NaN   NaN   NaN   2     2     NaN   NaN   NaN   NaN
    NaN   NaN   NaN     2   1     2     NaN   NaN   NaN   NaN   2     1     2     NaN   NaN   NaN
    NaN   NaN   2     1     1     2     NaN   NaN   NaN   NaN   2     1     1     2     NaN   NaN
    NaN   2     1     1     2     NaN   NaN   NaN   NaN   NaN   NaN   2     1     1     2     NaN
    2     1     1     1     2     2     2     2     2     2     2     2     1     1     1     2
    1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
    1     1     1     1     1     1     1     1     1     1     1     1     1     1     1     1
    2     1     1     1     2     2     2     2     2     2     2     2     1     1     1     2
    NaN   2     1     1     2     NaN   NaN   NaN   NaN   NaN   NaN   2     1     1     2     NaN
    NaN   NaN   2     1     1     2     NaN   NaN   NaN   NaN   2     1     1     2     NaN   NaN
    NaN   NaN   NaN   2     1     2     NaN   NaN   NaN   NaN   2     1     2     NaN   NaN   NaN
    NaN   NaN   NaN   NaN   2     2     NaN   NaN   NaN   NaN   2     2     NaN   NaN   NaN   NaN
    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
    ];

    function WinBtnMotionCallback0(src,evnt)
        cpos  = get(hax,'CurrentPoint');

        if abs(pos-cpos(1)) <= 3
            set(src,'pointer','custom','pointershapecdata',myptr,'pointershapehotspot',[9 9]);
        else
            set(src,'pointer','arrow');
        end
    end


    function WinBtnDownCallback(src,evnt)
        cpos  = get(hax,'CurrentPoint');
        if abs(pos-cpos(1)) <= 3
            set(src,'WindowButtonMotionFcn',@WinBtnMotionCallback)
            set(src,'WindowButtonUpFcn',@WinBtnUpCallback)
        end

        function WinBtnMotionCallback(src,evnt)
            cpos = get(hax,'CurrentPoint');
            cpos = min(max(1,cpos(1)),m(2));
            pos = round(cpos);
            set(hcurtain,'xdata',pos(1)*[1 1]);
            C(:,1:pos) = 1; C(:,pos+1:end)=0;
            set(hR,'alphadata',C);
            set(hT,'alphadata',1-C);
        end

        function WinBtnUpCallback(src,evnt)
            set(src,'WindowButtonMotionFcn',@WinBtnMotionCallback0)
            set(src,'WindowButtonUpFcn','')
        end
    end


end

