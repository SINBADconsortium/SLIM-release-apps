function [A n d o] = ReadAllData(file)
	%% [A n d o] = ReadAllData(file)
	% Read data in different type file
	% Right now, only support .rsf and .mat


[pathstr,name,ext] = fileparts(file);

switch ext
    case '.rsf' 
		[A n d o] = rsf_read_all(file);
		
	case '.mat'
		load(file)
		A = DataStructure.Data;
		n = DataStructure.n;
		o = DataStructure.o;
		d = DataStructure.d;
		
	otherwise
		error('Wrong data format')
end