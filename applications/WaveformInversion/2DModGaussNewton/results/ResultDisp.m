function ResultDisp(varargin)
% You can view FWI result with this function, top image is result, while 
% the bottom one is update.
%
% use: ResultDisp(name1,name2,...,videoRecord,'modeltype')
%
% Inputs:
%	name: name of result mat file
%	videoRecord: if videoRecord = 1, record a video; if 0, not
%	modeltype: modeltype = 'bg', BG compass model; 'camenbert', camenbert model
%
% Author: Xiang Li
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%         
% Date: 02, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.
%

if strcmp(varargin{end},'bg')
	vmin = 1480;
	vmax = 4600;
	varargin(end) = [];
elseif strcmp(varargin{end},'camenbert')
	vmin = 2250;
	vmax = 2750;
	varargin(end) = [];
else
	vmin = 1500;
	vmax = 3000;
end

Tint = .3;
fps  = 4;

if  varargin{end} == 1 | varargin{end} == 0
    video         = varargin{end};
    varargin(end) = [];
else
    video         = 0;
end


Numimgs  = length(varargin);
fig      = figure(1);

for m = 1:2*Numimgs
    subplot(2,Numimgs,m)
end


disp('Press any Key to continue')
pause

videoname = ['Result_of_'];
for m = 1:Numimgs
    name = varargin{m};
    if m == 1
        videoname = [videoname,varargin{m}];
    else
        videoname = [videoname,'and',varargin{m}];
    end
    eval(['data',num2str(m),'= load(name);'])
end

Movlength = size(data1.results,3);

if video
   aviobj = avifile([videoname,'.avi'],'fps',fps);
end


for m = 1:Movlength
    for n = 1:Numimgs
        subplot(2,Numimgs,n)
        eval(['imagesc(data',num2str(n),'.results(:,:,m))'])
        caxis([vmin vmax])
        colorbar;title([varargin{n}])
        subplot(2,Numimgs,Numimgs+n)
        eval(['imagesc(data',num2str(n),'.updates(:,:,m))'])
        colorbar;title(['update',num2str(m)])
    end
    pause(Tint)
    
    if video
        F = getframe(fig);
        aviobj = addframe(aviobj,F);
    end
    
end


if video
   aviobj = close(aviobj);
end
