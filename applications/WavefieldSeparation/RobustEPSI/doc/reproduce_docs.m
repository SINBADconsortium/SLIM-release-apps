% Run this to generate the html documentation from the Matlab source

curdir = pwd;

options.outputDir = curdir;
options.format = 'html';
options.stylesheet = 'slim.xsl';
options.createThumbnail = false;
options.showCode = true;

options.evalCode = false;
publish('index',options)

options.evalCode = true;
publish('example',options)

close all
open('index.html')