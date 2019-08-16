% testing window functions
clear all;
%N=[10 1 20]; % range of vector sizes
N=[10 1 1024]; % range of vector sizes
P=13;         % range of processors
H=13;         % range of half-overlap
T=13;        % tolerrance 1e-T

disp('start test');
fprintf('oplWindows1D: n=[%d:%d:%d] p=[2:%d] h=[0:%d] T=1e-%d\n',N(1),N(2),N(3),P,H,T);

for n=N(1):N(2):N(3)
    for p=2:P
        for h=0:H
    
            try
                [ m osf ysf xsf ] = pSPOT.pWindow.funWindow1DShape( n, p, h );
                %fprintf('\ttest: n=%d p=%d h=%d\n',n,p,h);
            catch pr
                %fprintf('invalid shape: n=%d p=%d h=%d\n',n,p,h);
                %disp(pr.message);
                continue;
            end

            try
                A=oplWindow1Dfd(n,p,h);
                try
                    A.utest(1);
                catch xy
                    fprintf('FD: n=%d m=%d p=%d h=%d\n',n,m,p,h);
                    disp(xy.message);
                    %disp(full(A))
                    %disp(full(A'))
                    %disp(full(A'*A))
                end
            catch AB
                fprintf('FD: n=%d m=%d p=%d h=%d\n',n,m,p,h);
                disp(AB.message);
            end
    
            try
                A=oplWindow1Dtpr(n,p,h);
                try
                    A.utest(1);
                catch xy
                    fprintf('TPR: n=%d m=%d p=%d h=%d\n',n,m,p,h);
                    disp(xy.message);
                    %disp(full(A))
                    %disp(full(A'))
                    %disp(full(A'*A))
                end
            catch AB
                fprintf('TPR: n=%d m=%d p=%d h=%d\n',n,m,p,h);
                disp(AB.message);
            end

            try
                A=oplWindow1Davg(n,p,h);
                try
                    A.utest(1);
                catch xy
                    fprintf('AVG: n=%d m=%d p=%d h=%d\n',n,m,p,h);
                    disp(xy.message);
                    %disp(full(A))
                    %disp(full(A'))
                    %disp(full(A'*A))
                end
            catch AB
                fprintf('AVG: n=%d m=%d p=%d h=%d\n',n,m,p,h);
                disp(AB.message);
            end

        end
    end
end

disp('end test');
