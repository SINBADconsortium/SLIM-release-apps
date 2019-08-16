function plot_group(A,x)
n = size(A,2);
linestyles = cellstr(char(':',':','-.','--',':',':','-.','--','-',':','-',':',...
'-.','--','-',':','-.','--','-',':','-.'));
 
MarkerEdgeColors=jet(n);  % n is the number of different items you have
Markers=['o','x','v','*','s','d','+','^','<','>','p','h','.',...
'+','*','o','x','^','<','h','.','>','p','s','d','v',...
'o','x','+','*','s','d','v','^','<','>','p','h','.'];
 
% [...]
 
for i=1:n
    if nargin < 2
        plot(A(:,i),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:),'linewidth',2);
    else
        plot(x,A(:,i),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:),'linewidth',2);
    end
    hold on
end