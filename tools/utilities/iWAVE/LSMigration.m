function [dm] = LSMigration(m,dU,model,opt)

% Time domain L1 least square migration code
% Note, right now we suggest only use data decompostion, not do the domain decomposition,
% because the domain decompostion of iwave++ seems not efficeiently

itermax  = 10;
L1M      = 1;
tol      = 1e-6;     



if isfield(opt,'itermax')
    itermax = opt.itermax;
end


if isfield(opt,'L1M')
    L1M     = opt.L1M;
end

if isfield(opt,'tol')
    tol     = opt.tol;
end



J        = oppDFd2DT(m,model,opt);

if L1M == 1
   opt.itermax = itermax;
   opt.tol     = tol; 
   C = opCurvelet(model.n(1),model.n(2));
   A = J*C';  
   dc = L1Solver(A,dU,C,opt);
   dm = C'*dc;
   
else
   if isfield(opt,'x0')
    x0 = opt.x0;
   else x0 = [];
   end
   %dm  = lsqrSOL(J,dU,opt);
   dm   = pcg(J'*J,J'*dU,tol,itermax,[],[],x0);
   %dm   = lsqr(J,dU,tol,itermax);
end   
