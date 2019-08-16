function A = opHelm2D9pt(vel,X,Y,nx_border,ny_border,f)
  % OPHELM2D9pt
  %
  % This function is used to construct matrix A from discretization of     
  % the 2D Helmholtz equation                                              
  %                                                                        
  %    (d^2/dx^2 + d^2/dy^2 + k^2) p = f                                   
  %                                                                        
  % using the 9-point finite difference scheme numerically optimized for minimal dispersion (Cho et. al., 97)
  % 
  % op = opHelm2D(vel,X,Y,nx_border,ny_border,f)
  % 
  % arguments:
  %
  %         vel         The velocity information as a matrix (will crash for 1D)
  %         X           Physical size of domain in X direction
  %         Y           Physical size of domain in Y direction
  %         nx_border   Number of damping boundry grids to use in x direction (typical 5~10% of domain)
  %         ny_border   Number of damping boundry grids to use in y direction (typical 5~10% of domain)
  %                     NOTE: total size of domain is (ny + 2*ny_border) by (nx + 2*nx_border)
  %         f           physical frequency of the H2 operator
  %
  % Author: Tim T.Y. Lin
  %
  %         AUG 18, 2009: Fixed Laplacian term coefficients due to typo on original paper
  

  % ----------------------- Initiation process

    % ah, ch, dh are parameters for the 9-pt stencil (see Jo et al), ah=1, ch=1, dh=0 degrades to defauls 5pt stencil
    ah = 0.5461;
    ch = 0.6248;
    dh = 0.09381;

% ah=1; 
% ch=1;
% dh=0; 


    [ny nx] = size(vel);
    x_length = X;
    y_length = Y;
    
    
    % rearrange velocity image into the subroutine's natural order:
    % 
    % ^
    % |     y-axis
    % |
    % |
    % ------------->  x-axis
    
    c = flipdim(vel,1);
    c = permute(c,[2,1]);
    
    fref = f;                 % reference wavelength at the top layer
    hx = x_length/(nx-1);           % mesh size in x-direction
    hy = y_length/(ny-1);           % mesh size in y-direction

    if (abs((hx - hy)/hx) > 1e-4)
        errstr = ['gridsize mismatch between x and y coordinates: hx/hy = ' num2str(abs(hx./hy))];
        error(errstr);
    end

    % enforce gridsize in y direction
    hy = hx;
    y_length = hx * (ny - 1);

    % global gridsize
    h = hx;
  
    % calculate border length
    if (ny_border == 0 || nx_border == 0)
        error('opHelm2D9pt needs a finite damping border size!')
    end
    
    if (ny_border < 1)   % length is a fraction
        ny_border = ceil(ny_border * ny);
    end
    
    ny_border_length = h * (ny_border + 1);
    
    if (nx_border < 1)
        nx_border = ceil(nx_border * nx);
    end
    
    nx_border_length = h * (nx_border + 1);
    
    % Declare some numbers
    nxy = nx * ny;
    w = 2*pi*f;
    a_fac = 1;
    
    ny_tot = ny + ny_border*2;
    nx_tot = nx + nx_border*2;
    
    % Laplacian term coefficients
    r1 = - (1 - ah) * (0.5 / h^2);
    r2 = - ah * (1 /h^2);
    r3 = - (1 - ah) * (0.5 /h^2);
    r4 = - ah * (1 / h^2);
    r5 = (ah * (4 / h^2)) + ((1 - ah) * (2 / h^2));
    r6 = - ah * (1 / h^2);
    r7 = - (1 - ah) * (0.5 / h^2);
    r8 = - ah * (1 / h^2);
    r9 = - (1 - ah) * (0.5 / h^2);
    
    % K2 term
    s1 = (1 - ch - 4*dh) / 4;
    s2 = dh;
    s3 = (1 - ch - 4*dh) / 4;
    s4 = dh;
    s5 = ch;
    s6 = dh;
    s7 = (1 - ch - 4*dh) / 4;
    s8 = dh;
    s9 = (1 - ch - 4*dh) / 4;
    
    % Second-order term in the boundry layers
    so2 = - 1 / h^2;
    so4 = - 1 / h^2;
    so5 = 4 / h^2;
    so6 = - 1 / h^2;
    so8 = - 1 / h^2;

    xlist = zeros(nxy*9,3);
    

    % Start grid definiton
    % 
    % This subroutine builds the finite difference matrix of Helholtz operator.
    % The matrix is stored in the 9-point stencil form, as follows:
    % 
    %        r7     r8     r9     s7    s8     s9          a7    a8    a9
    %            \   |  /             \  |  /                \   |  /
    %        r4 -   r5 -   r6  +  s4 -  s5  -  s6     =    a4 -  a5  - a6    
    %            /   |  \             /  |  \                /   |  \
    %        r1     r2     r3     s1    s2     s3          a1    a2    a3
    % 
    %      [_________________]   [---------------]        [---------------]
    %            Laplacian            k^2 term               Helmholtz  
    % 
    % 
    %  similarily I'm going to divide up the domain into 9 parts:
    %
    %
    %        d7     d8     d9  
    %            \   |  /      
    %        d4 -   d5 -   d6   
    %            /   |  \      
    %        d1     d2     d3  
    %                          
    %      [_________________] 
    %            Domain    
    %
    % The d5 part forms the entirety of the physical regieon. All the other parts are
    % the damping boundry layers, positioned as shown. They must be dealt separately
    % because the damping layer increases the complex velocity linearly in the outward
    % direction.

    
    % Firstly I divide up the index over the domain into its respective parts:
    % (remember, origin in this grid index is at the bottom left)
    d1_x_range = [1:nx_border];
    d1_y_range = [1:ny_border];
    
    [d1_y d1_x] = meshgrid(d1_y_range,d1_x_range);
    d1_x = d1_x(:);
    d1_y = d1_y(:);

    d2_x_range = [(nx_border + 1):(nx_border + nx)];
    d2_y_range = [1:ny_border];

    [d2_y d2_x] = meshgrid(d2_y_range,d2_x_range);
    d2_x = d2_x(:);
    d2_y = d2_y(:);

    d3_x_range = [(nx_border + nx + 1):nx_tot];
    d3_y_range = [1:ny_border];

    [d3_y d3_x] = meshgrid(d3_y_range,d3_x_range);
    d3_x = d3_x(:);
    d3_y = d3_y(:);

    d4_x_range = [1:nx_border];
    d4_y_range = [(ny_border + 1):(ny_border + ny)];

    [d4_y d4_x] = meshgrid(d4_y_range,d4_x_range);
    d4_x = d4_x(:);
    d4_y = d4_y(:);

    d5_x_range = [(nx_border + 1):(nx_border + nx)];
    d5_y_range = [(ny_border + 1):(ny_border + ny)];

    [d5_y d5_x] = meshgrid(d5_y_range,d5_x_range);
    d5_x = d5_x(:);
    d5_y = d5_y(:);

    d6_x_range = [(nx_border + nx + 1):nx_tot];
    d6_y_range = [(ny_border + 1):(ny_border + ny)];

    [d6_y d6_x] = meshgrid(d6_y_range,d6_x_range);
    d6_x = d6_x(:);
    d6_y = d6_y(:);

    d7_x_range = [1:nx_border];
    d7_y_range = [(ny_border + ny + 1):ny_tot];

    [d7_y d7_x] = meshgrid(d7_y_range,d7_x_range);
    d7_x = d7_x(:);
    d7_y = d7_y(:);
    
    d8_x_range = [(nx_border + 1):(nx_border + nx)];
    d8_y_range = [(ny_border + ny + 1):ny_tot];
    
    [d8_y d8_x] = meshgrid(d8_y_range,d8_x_range);
    d8_x = d8_x(:);
    d8_y = d8_y(:);

    d9_x_range = [(nx_border + nx + 1):nx_tot];
    d9_y_range = [(ny_border + ny + 1):ny_tot];
    
    [d9_y d9_x] = meshgrid(d9_y_range,d9_x_range);
    d9_x = d9_x(:);
    d9_y = d9_y(:);
    
    % Now make and inject coefficients into the stencil grid
    a = zeros(nx_tot*ny_tot,9);
        
    % d1
    d1_ind = sub2ind([nx_tot,ny_tot],d1_x,d1_y);
    v = w/c(1,1);
    alf = a_fac*(h*(nx_border-d1_x+1)/nx_border_length).^2 + a_fac*(h*(ny_border-d1_y+1)/ny_border_length).^2;
    a(d1_ind,1) = 0;
    a(d1_ind,2) = so2;
    a(d1_ind,3) = 0;
    a(d1_ind,4) = so4;
    a(d1_ind,5) = so5 - (1-i*alf).*(v.^2);
    a(d1_ind,6) = so6;
    a(d1_ind,7) = 0;
    a(d1_ind,8) = so8;
    a(d1_ind,9) = 0;
    
    % d2
    d2_ind = sub2ind([nx_tot,ny_tot],d2_x,d2_y);
    v = w./c(d2_x - nx_border,1);
    alf = a_fac*(h*(ny_border-d2_y+1)/ny_border_length).^2;
    a(d2_ind,1) = 0;
    a(d2_ind,2) = so2;
    a(d2_ind,3) = 0;
    a(d2_ind,4) = so4;
    a(d2_ind,5) = so5 - (1-i*alf).*(v.^2);
    a(d2_ind,6) = so6;
    a(d2_ind,7) = 0;
    a(d2_ind,8) = so8;
    a(d2_ind,9) = 0;
    
    % d3
    d3_ind = sub2ind([nx_tot,ny_tot],d3_x,d3_y);
    v = w./c(nx,1);
    alf = a_fac*(h*(d3_x-nx-nx_border)/nx_border_length).^2 + a_fac*(h*(ny_border-d3_y+1)/ny_border_length).^2;
    a(d3_ind,1) = 0;
    a(d3_ind,2) = so2;
    a(d3_ind,3) = 0;
    a(d3_ind,4) = so4;
    a(d3_ind,5) = so5 - (1-i*alf).*(v.^2);
    a(d3_ind,6) = so6;
    a(d3_ind,7) = 0;
    a(d3_ind,8) = so8;
    a(d3_ind,9) = 0;
    
    % d4
    d4_ind = sub2ind([nx_tot,ny_tot],d4_x,d4_y);
    v = w./c(1,d4_y-ny_border); v = v(:);
    alf = a_fac*(h*(nx_border-d4_x+1)/nx_border_length).^2;
    a(d4_ind,1) = 0;
    a(d4_ind,2) = so2;
    a(d4_ind,3) = 0;
    a(d4_ind,4) = so4;
    a(d4_ind,5) = so5 - (1-i*alf).*(v.^2);
    a(d4_ind,6) = so6;
    a(d4_ind,7) = 0;
    a(d4_ind,8) = so8;
    a(d4_ind,9) = 0;
    
    % d5
    d5_ind = sub2ind([nx_tot,ny_tot],d5_x,d5_y);
    d5_c_ind = sub2ind([nx,ny],d5_x-nx_border,d5_y-ny_border);
    v = w./c(d5_c_ind); v = v(:);
    alf = 0;
    a(d5_ind,1) = r1 - s1*(v.^2);
    a(d5_ind,2) = r2 - s2*(v.^2);
    a(d5_ind,3) = r3 - s3*(v.^2);
    a(d5_ind,4) = r4 - s4*(v.^2);
    a(d5_ind,5) = r5 - s5*(v.^2);
    a(d5_ind,6) = r6 - s6*(v.^2);
    a(d5_ind,7) = r7 - s7*(v.^2);
    a(d5_ind,8) = r8 - s8*(v.^2);
    a(d5_ind,9) = r9 - s9*(v.^2);
    
    % d6
    d6_ind = sub2ind([nx_tot,ny_tot],d6_x,d6_y);
    v = w./c(nx,d6_y-ny_border); v=v(:);
    alf = a_fac*((h*(d6_x-nx-nx_border))/nx_border_length).^2;
    a(d6_ind,1) = 0;
    a(d6_ind,2) = so2;
    a(d6_ind,3) = 0;
    a(d6_ind,4) = so4;
    a(d6_ind,5) = so5 - (1-i*alf).*(v.^2);
    a(d6_ind,6) = so6;
    a(d6_ind,7) = 0;
    a(d6_ind,8) = so8;
    a(d6_ind,9) = 0;
    
    % d7
    d7_ind = sub2ind([nx_tot,ny_tot],d7_x,d7_y);
    v = w./c(1,ny);
    alf = a_fac*(h*(nx_border-d7_x+1)/nx_border_length).^2 + a_fac*(h*(d7_y-ny-ny_border)/ny_border_length).^2;
    a(d7_ind,1) = 0;
    a(d7_ind,2) = so2;
    a(d7_ind,3) = 0;
    a(d7_ind,4) = so4;
    a(d7_ind,5) = so5 - (1-i*alf).*(v.^2);
    a(d7_ind,6) = so6;
    a(d7_ind,7) = 0;
    a(d7_ind,8) = so8;
    a(d7_ind,9) = 0;
    
    % d8
    d8_ind = sub2ind([nx_tot,ny_tot],d8_x,d8_y);
    v = w./c(d8_x-nx_border,ny); v=v(:);
    alf = a_fac*(h*(d8_y-ny-ny_border)/ny_border_length).^2;
    a(d8_ind,1) = 0;
    a(d8_ind,2) = so2;
    a(d8_ind,3) = 0;
    a(d8_ind,4) = so4;
    a(d8_ind,5) = so5 - (1-i*alf).*(v.^2);
    a(d8_ind,6) = so6;
    a(d8_ind,7) = 0;
    a(d8_ind,8) = so8;
    a(d8_ind,9) = 0;
    
    % d9
    d9_ind = sub2ind([nx_tot,ny_tot],d9_x,d9_y);
    v = w./c(nx,ny);
    alf = a_fac*(h*(d9_x-nx-nx_border)/nx_border_length).^2+ a_fac*(h*(d9_y-ny-ny_border)/ny_border_length).^2;
    a(d9_ind,1) = 0;
    a(d9_ind,2) = so2;
    a(d9_ind,3) = 0;
    a(d9_ind,4) = so4;
    a(d9_ind,5) = so5 - (1-i*alf).*(v.^2);
    a(d9_ind,6) = so6;
    a(d9_ind,7) = 0;
    a(d9_ind,8) = so8;
    a(d9_ind,9) = 0;
    
    
    % Boundary condition at edges  
    % 
    %   N, S, W, E : all radiation boundary conditions
    % 
    % N : 
    ix = 1:nx_tot;
    iy = ny_tot*ones(1,length(ix));
    ind = sub2ind([nx_tot,ny_tot],ix,iy);
    ind = ind(:);
    
    v_bord1 = (w/c(1,ny)) * ones(1,nx_border);
    v_mid   = w./c([1:nx],ny);
    v_bord2 = (w/c(nx,ny)) * ones(1,nx_border);
    v = [v_bord1(:); v_mid(:); v_bord2(:)];
    
    a(ind,5) = a(ind,5) - (1 ./ (1 + i.*v.*h))./(h^2);
    a(ind,7) = 0;
    a(ind,8) = 0;
    a(ind,9) = 0;
    
    % S : 
    ix = 1:nx_tot;
    iy = 1*ones(1,length(ix));
    ind = sub2ind([nx_tot,ny_tot],ix,iy);
    ind = ind(:);
    
    v_bord1 = (w/c(1,1)) * ones(1,nx_border);
    v_mid   = w./c([1:nx],1);
    v_bord2 = (w/c(nx,1)) * ones(1,nx_border);
    v = [v_bord1(:); v_mid(:); v_bord2(:)];
    
    a(ind,5) = a(ind,5) - (1 ./ (1 + i.*v.*h))./(h^2);
    a(ind,1) = 0;
    a(ind,2) = 0;
    a(ind,3) = 0;

    % W : 
    iy = 1:ny_tot;
    ix = 1*ones(1,length(iy));
    ind = sub2ind([nx_tot,ny_tot],ix,iy);
    ind = ind(:);
    
    v_bord1 = (w/c(1,1)) * ones(1,ny_border);
    v_mid   = w./c(1,[1:ny]);
    v_bord2 = (w/c(1,ny)) * ones(1,ny_border);
    v = [v_bord1(:); v_mid(:); v_bord2(:)];
    
    a(ind,5) = a(ind,5) - (1 ./ (1 + i.*v.*h))./(h^2);
    a(ind,1) = 0;
    a(ind,4) = 0;
    a(ind,7) = 0;

    % E :     
    iy = 1:ny_tot;
    ix = nx_tot*ones(1,length(iy));
    ind = sub2ind([nx_tot,ny_tot],ix,iy);
    ind = ind(:);
    
    v_bord1 = (w/c(nx,1)) * ones(1,ny_border);
    v_mid   = w./c(nx,[1:ny]);
    v_bord2 = (w/c(nx,ny)) * ones(1,ny_border);
    v = [v_bord1(:); v_mid(:); v_bord2(:)];
    
    a(ind,5) = a(ind,5) - (1 ./ (1 + i.*v.*h))./(h^2);
    a(ind,3) = 0;
    a(ind,6) = 0;
    a(ind,9) = 0;

    % return matrix to the original orientation of the velocity image
    a = reshape(a,[nx_tot,ny_tot,9]);
    a = permute(a,[2,1,3]);
    a = flipdim(a,1);
    a = reshape(a,[ny_tot,nx_tot,3,3]);
    a = permute(a,[1,2,4,3]);
    a = flipdim(a,3);
    a = reshape(a,[ny_tot,nx_tot,9]);

    %  NOTE: after the re-ordering, 'a' stencil is now in the orientation
    %  compatible with the MATLAB (:) operator, namely
    %
    %
    %        a1     a4     a7  
    %            \   |  /      
    %        a2 -   a5 -   a8   
    %            /   |  \      
    %        a3     a6     a9  
    %                          
    %      [_________________] 
    %            Domain
    
    a = reshape(a,[nx_tot*ny_tot,9]);
    sp_rowind = [];
    sp_colind = [];
    sp_vals = [];
    
    % == inject the stencil into sparse matrix format ==
    % build the coordinate/value vectors for the sparse matrix
    
    % main diagonal a5
    ind = find(a(:,5));
    sp_rowind = [sp_rowind; ind];
    sp_colind = [sp_colind; ind];
    sp_vals = [sp_vals; a(ind,5)];
    
    % a1
    ind = find(a(:,1));
    sp_rowind = [sp_rowind; ind];
    sp_colind = [sp_colind; ind-1-ny_tot];
    sp_vals = [sp_vals; a(ind,1)];
    
    % a2
    ind = find(a(:,2));
    sp_rowind = [sp_rowind; ind];
    sp_colind = [sp_colind; ind-ny_tot];
    sp_vals = [sp_vals; a(ind,2)];
    
    % a3
    ind = find(a(:,3));
    sp_rowind = [sp_rowind; ind];
    sp_colind = [sp_colind; ind+1-ny_tot];
    sp_vals = [sp_vals; a(ind,3)];
    
    % a4
    ind = find(a(:,4));
    sp_rowind = [sp_rowind; ind];
    sp_colind = [sp_colind; ind-1];
    sp_vals = [sp_vals; a(ind,4)];
    
    % a6
    ind = find(a(:,6));
    sp_rowind = [sp_rowind; ind];
    sp_colind = [sp_colind; ind+1];
    sp_vals = [sp_vals; a(ind,6)];
    
    % a7
    ind = find(a(:,7));
    sp_rowind = [sp_rowind; ind];
    sp_colind = [sp_colind; ind-1+ny_tot];
    sp_vals = [sp_vals; a(ind,7)];
    
    % a8
    ind = find(a(:,8));
    sp_rowind = [sp_rowind; ind];
    sp_colind = [sp_colind; ind+ny_tot];
    sp_vals = [sp_vals; a(ind,8)];
    
    
    % a9
    ind = find(a(:,9));
    sp_rowind = [sp_rowind; ind];
    sp_colind = [sp_colind; ind+1+ny_tot];
    sp_vals = [sp_vals; a(ind,9)];
    
    % keyboard
    
    A = sparse(sp_rowind, sp_colind, sp_vals);
    op = @(x,mode) opHelm2D9pt_intrnl(A,mode,x);

        
    function y = opHelm2D9pt_intrnl(A,mode,x)
    
        switch mode  
            case 0
                c =~isreal(A);
                [m n] = size(A);
                y = {m,n,[c,1,c,1],opinfo};
    
            case 1
                y = A * x;
    
            case 2
                y = A' * x;
        end
    end
 
end
        

