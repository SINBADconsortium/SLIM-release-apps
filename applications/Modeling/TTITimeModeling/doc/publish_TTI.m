clear all;
close all;
clc;

options.evalCode = true;
options.stylesheet = 'slim.xsl';
options.format = 'html';
options.outputDir = './';

file =publish('../scripts/Modelling_TTI_2D.m',options);
web(file);

file =publish('../scripts/Modelling_3D.m',options);
web(file);