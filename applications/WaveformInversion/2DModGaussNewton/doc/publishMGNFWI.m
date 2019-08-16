%% This script generate documentation for this software release
% 
options.evalCode = false;
options.stylesheet = 'slim.xsl';
options.format = 'html';
options.outputDir = './html/';


% camenbert example

publish('MGNFWI_camenbert.m',options);

% BG compress example

publish('MGNFWI_BG.m',options);


options.evalCode = true;

publish('examples.m',options);

publish('index.m',options);