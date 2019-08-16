% current directory
curdir = [pwd '/'];
basedir = curdir(1:end-length('doc/'));
scriptdir = [basedir 'scripts/'];
datadir = [basedir 'data/'];
resultsdir = [basedir 'results/'];

% initialize empty options struct to delete any previous values
options = [];

% output directory
options.outputDir = './HTML';

% location of stylesheet
options.stylesheet = [curdir 'slim.xsl'];

% format
options.format      = 'html';

% no thumbnails
options.createThumbnail = false;

% whether to show that code
options.showCode = true;


% publish files without running the code
options.evalCode = false;
publish('index.m',options);
publish([scriptdir 'camembert.m'],options);
publish([scriptdir 'bg2_batch.m'],options);


% Publish file and run the code, producing figures
options.evalCode = true;

publish('results_2drobustbatching.m',options);

% close figures
close all;

