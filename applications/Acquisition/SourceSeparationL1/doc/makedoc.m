% Current directory
curdir = pwd;

% Initialize empty options struct to delete any previous values
options = [];

% Output directory
options.outputDir = curdir;

% Location of stylesheet
options.stylesheet = [curdir '/slim.xsl'];

% Format
options.format = 'html';

% No thumbnails
options.createThumbnail = false;

% Show the code
options.showCode = true;

% Publish some files without running the code
options.evalCode = false;
publish('index.m',options);

% Publish a file and run the code
options.evalCode = true;
publish('example.m',options);

% Close figures
close all;

