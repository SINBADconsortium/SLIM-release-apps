%% This script generate documentation for this software release
% 
options.evalCode = false;
options.stylesheet = 'slim.xsl';
options.format = 'html';
options.outputDir = './html/';




options.evalCode = true;

publish('examples.m',options);

publish('index.m',options);