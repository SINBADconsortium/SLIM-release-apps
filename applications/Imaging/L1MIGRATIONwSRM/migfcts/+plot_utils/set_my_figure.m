function set_my_figure(ah,fontsize)
% ah is current axis handle, use gca to get it
% by Ning Tu, SLIM, 2011
if nargin < 1
    ah = gca;
    fontsize = 12;
elseif nargin < 2
    fontsize = 12;
end
set(ah,'FontName','Arial','FontSize',fontsize,'FontWeight','normal','TickDir','out')
set(get(ah,'XLabel'),'FontName','Arial','FontSize',fontsize,'FontWeight','normal')
set(get(ah,'YLabel'),'FontName','Arial','FontSize',fontsize,'FontWeight','normal')
set(get(ah,'Title'),'FontName','Arial','FontSize',fontsize,'FontWeight','normal')