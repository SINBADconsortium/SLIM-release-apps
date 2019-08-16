% To reproduce the documentation
option.format = 'html';
option.stylesheet = 'slim.xsl';
% no thumbnails
options.createThumbnail = false;
% whether to show that code
options.showCode = true;
options.evalCode = true;

publish('index',option);

options.evalCode = true;
publish('interp4Dview',option);

close all;