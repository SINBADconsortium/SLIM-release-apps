function test_geomCG
   
   if exist('geomCG.m','file') == 0
       warning('geomCG must be installed in order to run this test, skipping');
   else       
       n = 100; d = 3; r = 20*ones(1,d); p = 0.1;
       dims = n*ones(1,d);
       
       % Check correctness + performance 
       % We pass the same data to both geomCG and fitTucker and let them run for the same # of iterations
       % The generated tensor should be the same for both methods
       
       A = makeRandTensor( dims, r );
       % Generate data
       subs = makeOmegaSet( dims, round(p*prod(dims)) );
       vals = getValsAtIndex(A, subs);

       A_Omega = sptensor( subs, vals, dims, 0);                    
       
       % random initial guess:
       X_init = makeRandTensor( dims, r );
       
       subs_1d = sub2ind(dims,subs(:,1),subs(:,2),subs(:,3));
       e = false(prod(dims),1);
       vals_full = zeros(prod(dims),1);
       e(subs_1d) = true;
       vals_full(subs_1d) = vals;
       
       % geomCG optimization params
       maxIter = 10;
       opts = struct( 'maxiter', maxIter, 'tol', 1e-9 );
       x0 = [];
       for i=1:d
           x0 = [x0; vec(X_init.U{i})]; 
       end
       x0 = [x0; vec(X_init.core.data)];
       
       tic;
       Xgeom = geomCG( A_Omega, X_init, [], opts);
       geomcgTime = toc;       
       
       % fitTucker optimization
       tic;
       [Udense,Bdense] = fitTucker( e, vals_full, dims, r ,'x0',x0,'maxiter',maxIter-1 , 'maxLS',1 );
       Xdense = ttm(Bdense,Udense);
       denseTime = toc;       
       
       % Results
       % Roundoff errors change the result the more iterations are performed
       err = norm(vec(double(Xgeom)) - vec(Xdense))/norm(vec(Xdense));
       assertTrue(err < 1e-2 );
       disp(['geomCG time : ' num2str(geomcgTime) 's']); 
       disp(['dense Tucker time : ' num2str(denseTime) 's']);
   end