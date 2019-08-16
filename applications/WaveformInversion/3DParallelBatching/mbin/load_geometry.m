function  model  = load_geometry( model,model_name,numsrc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nx = model.n(1); ny = model.n(2); nz = model.n(3);
[x,y,z] = odn2grid(model.o,model.d,model.n);
switch model_name
  case 'edam'
    model.xsrc = x(round(linspace(3,nx-3,numsrc)));
    model.ysrc = y(round(linspace(3,ny-3,numsrc)));
    model.zsrc = min(z)+model.d(3);
    model.xrec = x(3:1:end-2);
    model.yrec = y(3:1:end-2);
    model.zrec = max(z)-model.d(3);
    model.t0 = 0;
    model.f0 = 10;
  
  case {'overthrust','overthrust_small'}
    model.xsrc = x(round(linspace(3,nx-3,numsrc)));
    model.ysrc = y(round(linspace(3,ny-3,numsrc)));
    model.zsrc = min(z)+model.d(3);
    model.xrec = x(3:1:end-2);
    model.yrec = y(3:1:end-2);
    model.zrec = min(z)+model.d(3);
    
    model.t0   = 0;                
    model.f0 = 10;
end
end

