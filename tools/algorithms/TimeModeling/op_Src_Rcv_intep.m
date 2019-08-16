classdef op_Src_Rcv_intep < opSpot
	%% op = op_Src_Rcv_intep(mmz,mmx,mmy,z,x,y)
	% This function helps you finds values of some points in a 2D/3D data cube
	% Those points can be not on the grid.
	% 
	%	Input:
	%		mmz: grid coordinate of 2D/3D cube along z, length(mmz) = size(data,1);
	%		mmx: grid coordinate of 2D/3D cube along x, length(mmx) = size(data,2);
	%		mmy: grid coordinate of 2D/3D cube along y, length(mmy) = size(data,3);
	%			NOTICE, Grid needs to have same space interval along each dimension.
	%		[z,x,y]: n x 3 matrix, each row is the coordinate of the point where
	%				you want the value.
	%
	%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = protected)
   		funHandle = []; 
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = op_Src_Rcv_intep(mmz,mmx,mmy,z,x,y)
		   
		   if nargin < 6, y = ones(size(x));end
		   z-x-y;% for check, z, x, y need to be same length
		   np = length(z);
		   
		   siz = [length(mmz),length(mmx),length(mmy)];
		
		   
		   if length(mmy) == 1

	 		   A = spalloc(np,prod(siz),np*16);
		   
				for mm = 1:np
					
					x_id  = length(find(x(mm)-mmx>=0));
					z_id  = length(find(z(mm)-mmz>=0));
					x_id_grid = [floor(x_id)-1,floor(x_id),floor(x_id)+1,floor(x_id)+2];
					z_id_grid = [floor(z_id)-1,floor(z_id),floor(z_id)+1,floor(z_id)+2];
					
					ceos_x = lagrange_ceo(x(mm),mmx(x_id_grid));
					ceos_z = lagrange_ceo(z(mm),mmz(z_id_grid));
					ceos   = ceos_z*ceos_x';
					
					[XX,ZZ]   = meshgrid(x_id_grid,z_id_grid);
					ind		  = sub2ind(siz,ZZ,XX);
					

					A(mm,ind(:))  = ceos(:);
				end	
				
				
		   else
	
	 		   A = spalloc(np,prod(siz),np*64);
		   
				for mm = 1:np
					
					x_id  = length(find(x(mm)-mmx>=0));
					z_id  = length(find(z(mm)-mmz>=0));
					y_id  = length(find(y(mm)-mmy>=0));
					x_id_grid = [floor(x_id)-1,floor(x_id),floor(x_id)+1,floor(x_id)+2];
					z_id_grid = [floor(z_id)-1,floor(z_id),floor(z_id)+1,floor(z_id)+2];
					y_id_grid = [floor(y_id)-1,floor(y_id),floor(y_id)+1,floor(y_id)+2];
					
					ceos_x = lagrange_ceo(x(mm),mmx(x_id_grid));
					ceos_z = lagrange_ceo(z(mm),mmz(z_id_grid));
					ceos_y = lagrange_ceo(y(mm),mmy(y_id_grid));
					ceos   = kron(ceos_y',kron(ceos_x',ceos_z'));
			
					[XX,ZZ,YY]   = meshgrid(x_id_grid,z_id_grid,y_id_grid);
					ind		  = sub2ind(siz,ZZ,XX,YY);
					

					A(mm,ind(:)) = ceos(:);
				end	
			
		   end
		   

		   % Construct operator
         
  			fun = @(x,mode) op_Src_Rcv_intep_intrnl(siz,A,x,mode);
            % Construct operator
            op = op@opSpot('Src_Rcv_intep', np, prod(siz));
            op.cflag     = 1;
            op.funHandle = fun;
       end % Constructor
 
    end % Methods
       
 
    methods ( Access = protected )
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    % Multiply
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    function y = multiply(op,x,mode)
	       y = op.funHandle(x,mode);
	    end % Multiply          

	 end % Methods
   
end % Classdef



function y = op_Src_Rcv_intep_intrnl(siz,A,x,mode)
	if (mode == 1)
		% if size(data_index,1) == 1
		% 	y = x(data_index)'.* data_weight;
		% else
		% 	y = x(data_index).* data_weight;
		% end
		% y = sum(y,2);
		y = A * x;
		
	elseif (mode == 2)  
		% y = zeros(siz);
		% for mm = 1:size(data_weight,1);
		% 	y(data_index(mm,:)) = y(data_index(mm,:)) + repmat(x(mm),1,size(data_weight,2)).*data_weight(mm,:);
		% end
		% y = y(:);
		y = A' * x;
		
	end
end