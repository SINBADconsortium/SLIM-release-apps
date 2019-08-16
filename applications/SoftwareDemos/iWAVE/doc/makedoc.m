% current directory
curdir = pwd;

% initialize empty options struct to delete any previous values
options = [];

% output directory
options.outputDir = './HTML';

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
publish('makedoc.m',options);

% publish a file and run the code
options.evalCode = true;

publish('example.m',options);

% close figures
close all;
% open file for viewing
open('./HTML/index.html');
open('./HTML/example.html');
open('./HTML/makedoc.html');
