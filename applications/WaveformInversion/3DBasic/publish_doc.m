doc_dir = [pwd '/doc/'];
out_dir = [pwd '/doc/html/'];

options.stylesheet = 'slim.xsl';
options.format = 'html';
options.outputDir = out_dir;
options.evalCode = true;

publish([doc_dir 'example_edam1.m'],options);

options.evalCode = false;

publish([doc_dir 'index.m'],options);