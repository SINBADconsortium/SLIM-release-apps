% current directory
curdir = pwd;

% initialize empty options struct to delete any previous values
options = [];

% output directory
options.outputDir = curdir;

% location of stylesheet
options.stylesheet = [curdir '/slim.xsl'];

% format
options.format      = 'html';

% no thumbnails
options.createThumbnail = false;

% whether to show that code
options.showCode = true;


% publish some files without running the code
options.evalCode = false;

publish('index.m',options);
publish('wL1_freq_MH.m',options);
publish('wL1min.m',options);

% publish a file and run the code
options.evalCode = true;

publish('results.m',options);

% close figures
close all;
% % open file for viewing
% open('index.html')