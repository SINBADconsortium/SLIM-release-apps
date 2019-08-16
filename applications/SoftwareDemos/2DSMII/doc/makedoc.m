
opts = [];
opts.outputDir  = pwd;
opts.stylesheet = 'slim.xsl';
opts.format     = 'html';

%% scripts
opts.evalCode   = false;
opts.showCode   = true;

publish('../scripts/testing_adjoint.m',opts);
publish('../scripts/testing_analytic.m',opts);
publish('../scripts/testing_jacobian.m',opts);
publish('../scripts/testing_parallel.m',opts);
publish('../scripts/example_fwi.m',opts);
publish('../scripts/example_rtm.m',opts);
publish('../scripts/example_tomo.m',opts);

%% pages 
opts.evalCode   = true;
opts.showCode   = true;

publish('modeling',opts);
publish('testing',opts);
publish('examples',opts);
html_file = publish('index',opts);

%% show
close all;
open(html_file);