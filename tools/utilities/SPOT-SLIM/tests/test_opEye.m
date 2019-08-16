function test_suite = test_opEye
%test_opDiag  Unit tests for the opDiag operator
test_suite=buildFunctionHandleTestSuite(localfunctions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seed = setup
   seed = rng('default');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_opEye_mv(seed)
    nrm=[];
    for m=3:4
        for n=3:4
            for i=1:2
                %fprintf('m=%d n=%d i=%d\n',m,n,i);
                A=opEye(m,n);
                nrm=[nrm;dottest(A)];
            
                % forward
                b=rand(n,i);
                y=double(A*b);
                if m>=n
                    bb=[b; zeros(m-n,i)];
                else
                    bb=b(1:m,:);
                end
                nrm=[nrm;norm(y-bb)];
            
                % inverse
                b=rand(m,i);
                y=double(A'*b);
                if n>=m
                    bb=[b; zeros(n-m,i)];
                else
                    bb=b(1:n,:);
                end
                nrm=[nrm;norm(y-bb)];
            
            end
        end
    end
    assertEqual(norm(nrm),0)
end
