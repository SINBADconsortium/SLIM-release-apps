This folder include all the results which all saved as MAT format

View result:

you can use 'ResultDisp' to view the result interactively, or you can use following
	way.

1. go to subfolder, e.g 'example_BG'
2. load 'FWIresult*.mat' file into matlab. If you haven't run any 
		experiment, just run 'scons' for the preview result.
3. After load the data, u will get 3 variables, 'results', 'updates', 'resultinfo'.
		'results' and 'updates' are 3D data cube; first two axis is depth and lateral,
		while the 3rd axis is GN iteration. Then you can use matlab display functions
		like 'imagesc' to display 2D slides of the data cube. 
4. 'resultinfo' is a data container which contains convergency information at each
		iteration.
