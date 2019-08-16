%% This script generate documentation for this software release
% 
options.evalCode = false;
options.stylesheet = 'slim.xsl';
options.format = 'html';
options.outputDir = './html/';


% camenbert example


% % BG compress example
%
publish('demo2D_gradient_test',options);


options.evalCode = true;

publish('TimeModeling_example.m',options);

publish('index.m',options);