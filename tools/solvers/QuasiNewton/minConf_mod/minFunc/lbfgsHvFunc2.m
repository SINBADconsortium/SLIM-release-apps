function Hv = lbfgsHvFunc2(v,Hdiag,N,M)
%size(M)
%cond(M)
if cond(M)>(1/(eps*1e3))
    pr = ssbin(M,500);
    L=spdiags(pr,0,length(pr),length(pr));
    Hv = v/Hdiag - N*L*((L*M*L)\(L*(N'*v)));
    
else
    Hv = v/Hdiag - N*(M\(N'*v));
end