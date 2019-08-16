function WriteAllData(file,Data,n,d,o)
	%%  WriteAllData(file,Data,n,d,o)
	% Wrte data into different type file
	% Right now, only support .rsf and .mat
	 
	 
	[pathstr,name,ext] = fileparts(file);

	switch ext
	    case '.rsf' 
			rsf_write_all(file,{'out=stdout'},Data,d,o);
		
		case '.mat'
			DataStructure.Data = Data;
			DataStructure.n    = n;
			DataStructure.o    = o;
			DataStructure.d    = d;
			save(file,'DataStructure','-v7.3');
		
		otherwise
			error('Wrong data format')
	end
