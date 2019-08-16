function [f]=extender2(A,transp_flag);

M=length(A(:,1));
N=length(A(1,:));

if (nargin >1) & strcmp (transp_flag,'transp')
    f=A(51:end-50,51:end-50);
else
    A=[repmat(A(:,1),1,50) A repmat(A(:,end),1,50)];
    f=[repmat(A(1,:),50,1); A; repmat(A(end,:),50,1)];
end
