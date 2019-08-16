function A = getLA(x1,x2)

%% Linear interpolation for source and receiver position
% find interior points
% Author : Mathias Louboutin
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earth, Ocean, and Atmosperic Sciences
%         The University of British Columbia
% Date : July, 2015
ik = find(x2 >= x1(1) & x2 <= x1(end));
% sizes
n1=length(x1);n2=length(x2);nk=length(ik);
% check
if ~nk
    A = [];
    return;
end

% initialize stuff
I = zeros(4*nk,1); J = I; S = I;
a=1;b=2;c=3;d=4;
l = 1;

% loop
for i = 1:nk
    k = ik(i);
    if x2(k)<x1(b)
        while (x2(k)<x1(b))&&(b-1>1)
            b=b-1;
        end
        a=b-1;c=b+1;d=c+1;
    elseif x2(k)>x1(c)
        while (x2(k)>x1(c))&&(c+1<n1)
            c=c+1;
        end
        a=c-2;b=c-1;d=c+1;
    end
    
    I(l)   = k; J(l)   = a; S(l)   = ((x2(k)-x1(b))*(x2(k)-x1(c))*(x2(k)-x1(d)))/((x1(a)-x1(b))*(x1(a)-x1(c))*(x1(a)-x1(d)));
    I(l+1) = k; J(l+1) = b; S(l+1) = ((x2(k)-x1(a))*(x2(k)-x1(c))*(x2(k)-x1(d)))/((x1(b)-x1(a))*(x1(b)-x1(c))*(x1(b)-x1(d)));
    I(l+2) = k; J(l+2) = c; S(l+2) = ((x2(k)-x1(b))*(x2(k)-x1(a))*(x2(k)-x1(d)))/((x1(c)-x1(b))*(x1(c)-x1(a))*(x1(c)-x1(d)));
    I(l+3) = k; J(l+3) = d; S(l+3) = ((x2(k)-x1(b))*(x2(k)-x1(c))*(x2(k)-x1(a)))/((x1(d)-x1(b))*(x1(d)-x1(c))*(x1(d)-x1(a)));
    l = l + 4;
end
% construct sparse matrix
A = sparse(I(1:l-1),J(1:l-1),S(1:l-1),n2,n1);
end