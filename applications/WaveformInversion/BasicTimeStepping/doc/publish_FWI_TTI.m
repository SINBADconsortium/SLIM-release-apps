close all;
clc;

options.evalCode = true;
options.stylesheet = 'slim.xsl';
options.format = 'html';
options.outputDir = './';

file =publish('../scripts/FWI_TTI_2D.m',options);
web(file);
