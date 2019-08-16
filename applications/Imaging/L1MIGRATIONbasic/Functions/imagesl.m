function h = imagesl(varargin)
%IMAGESC Scale data and display as image.
%   IMAGESL(...) is the same as IMAGE(...) except the data is scaled
%   to use the gray colormap.
%   
%   IMAGESC(...,CLIM) where CLIM = [CLOW CHIGH] can specify the
%   scaling.


clim = [];
switch (nargin),
  case 0,
    hh = image('CDataMapping','scaled');
  case 1,
    hh = image(varargin{1},'CDataMapping','scaled');
  case 3,
    xidex = max(size(varargin{1}));
    yidex = max(size(varargin{2}));
    varargin{3} = reshape(varargin{3},yidex,xidex);
    hh = image(varargin{:},'CDataMapping','scaled');
  otherwise,

    % Determine if last input is clim
    if isequal(size(varargin{end}),[1 2])
      str = false(length(varargin),1);
      for n=1:length(varargin)
        str(n) = ischar(varargin{n});
      end
      str = find(str);
      if isempty(str) || (rem(length(varargin)-min(str),2)==0),
        clim = varargin{end};
        varargin(end) = []; % Remove last cell
      else
        clim = [];
      end
    else
      clim = [];
    end
    hh = image(varargin{:},'CDataMapping','scaled');
end

% Get the parent Axes of the image
cax = ancestor(hh,'axes');

if ~isempty(clim),
  set(cax,'CLim',clim)
elseif ~ishold(cax),
  set(cax,'CLimMode','auto')
end

if nargout > 0
    h = hh;
end
colorbar;
colormap gray;
set(gca,'fontsize',18)
