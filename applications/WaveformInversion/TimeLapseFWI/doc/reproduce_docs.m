% To reproduce the documentation

% current directory
curdir = pwd;

% initialize empty options struct to delete any previous values
options = [];

% output directory
options.outputDir = curdir;

% location of stylesheet
options.stylesheet = 'slim.xsl';

% format
options.format = 'html';

% no thumbnails
options.createThumbnail = false;

% whether to show that code
options.showCode = true;

% publish some files without running the code
options.evalCode = false;
publish('index.m',options);

% publish files and run the code
options.evalCode = true;
publish('example_timelapse_BG.m',options);

% close figures
close all;

