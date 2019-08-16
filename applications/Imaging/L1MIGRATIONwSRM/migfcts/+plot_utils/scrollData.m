% scrollData(R,...);
function scrollData(R,varargin)
im      = { R , varargin{:} };
nimages = length(im);

m = size(R);
for k=1:nimages
    R = im{k};
    if size(R,3) ~= m(3),
        error('All image must have same 3rd dim');
    end;
    % true-color images
    if size(R,4) > 1 && (min(R(:))<0 || max(R(:))>1), 
        R=R-min(R(:));
        R=R/max(R(:));
        im{k} = R;
    end
end


imname  = cell(nimages,1);
for k=1:nimages
    if isempty(inputname(k)), imname{k} = sprintf('ARG %d',k);
    else                      imname{k} = sprintf('%s',inputname(k));
    end;
end


slice = round(m(3)/2);
p = round(sqrt(nimages));
q = ceil(nimages/p);
if (p-1)*q >= nimages, p=p-1; end


fig = figure;
clf;
set(fig, ...
    'WindowScrollWheelFcn',@figMouseWheelScroll, ...
    'KeyPressFcn',@figUpDownKeyScroll);

him  = zeros(nimages,1);
htxt = zeros(nimages,1);
for k =1:nimages,
    R = im{k};
    subplot(p,q,k);
    him(k)  = imagesc(squeeze(R(:,:,slice,:)));
    htxt(k) = title(sprintf('%d / %d x %d x %d   (%s)',slice,size(R,1),size(R,2),size(R,3),imname{k}));
end

if nimages == 1
    htext=text(1,1,'','fontsize',16,'color','y','verticalalignment','bottom');
    pos = get(gca,'CurrentPoint');
    row = round(pos(1,2));
    col = round(pos(1,1));
    set(fig,'WindowButtonMotionFcn',@WinBtnMotionCallback);
else
    set(fig,'WindowButtonMotionFcn',@sliderCallback);
end

hslider=uicontrol( ...
    'Style','slider', ...
    'Units','normalized', ...
    'position',[0.8 0 0.2 1 ], ...
    'value',slice, ...
    'min',1, ...
    'max',m(3));
set(hslider,'callback',@sliderCallback);
drawnow
    function sliderCallback(src,event)
        slice = round(get(hslider,'value'));
        drawSlice;
    end

    function figMouseWheelScroll(src,evnt)
        if evnt.VerticalScrollCount > 0
            if slice<m(3), slice = slice+1; drawSliceAndUpdateSlider; end
        elseif evnt.VerticalScrollCount < 0
            if slice>1,  slice = slice-1;  drawSliceAndUpdateSlider;  end
        end
    end

    function figUpDownKeyScroll(src,evnt)
        if strcmp(evnt.Key,'uparrow')
            if slice<m(3), slice = slice+1; drawSliceAndUpdateSlider; end
        elseif strcmp(evnt.Key,'downarrow')
            if slice>1,  slice = slice-1;   drawSliceAndUpdateSlider;  end
        end
    end

    function WinBtnMotionCallback(src,event)
        pos   = get(gca,'CurrentPoint');
        row = round(pos(1,2));
        col = round(pos(1,1));
        printtext
        sliderCallback(src,event);
    end


    function drawSliceAndUpdateSlider
        set(hslider,'value',slice);
        drawSlice;
    end

    function drawSlice
        for k =1:nimages,
            R = im{k};
            set(him(k),'cdata',squeeze(R(:,:,slice,:)));
            set(htxt(k),'string',sprintf('%d / %d x %d x %d   (%s)',slice,size(R,1),size(R,2),size(R,3),imname{k}));
        end
        printtext
        drawnow;
    end

    function printtext
        if nimages==1,
            if row>=1 && row<=m(1) && col>=1 && col<=m(2)
                str = sprintf('(%d,%d,%d) %s',row,col,slice,sprintf('%g ',R(row,col,slice,:)));
                set(htext,'pos',pos(1,[1 2]),'string',str);
                set(fig,'pointer','cross');
            else
                set(htext,'string','');
                set(fig,'pointer','arrow');
            end
        end
    end



end % scroll_wheel