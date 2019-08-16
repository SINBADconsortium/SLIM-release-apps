function [ m os ys xs ] = funWindow1DShape( n, p, h )
%funWindow1DShape is a support function for forward oplWindow1D* operators.
%
%   [ M OS YS XS ] = funWindow1DShape( N, P, H )
%
%   INPUT:
%      N = length of the input vector
%      P = number of processors
%      H = half of the overlap's size
%   OUTPUT:
%      M = length of the output vector
%      OS = (p,3) vector holding start, size, end indecies
%           of the default distribution of the input vector in every window
%      YS = (p,3) vector holding start, size, end indecies
%           of the default distribution of the output vector in every window
%      XS = (p,3) vector holding start, size, end indecies
%           of the input vector that will end up in every ouput window
%   NOTE:
%      This function assumes that the both the input and output vector
%      will follow the default distribution

    % check # of processors
    assert(p>1,'funWindow1DShape: number of processors has to be bigger then %d',p);

    % initialize shape vecs
    oo=zeros(p,1); od=zeros(p,1); oe=zeros(p,1);
    yo=zeros(p,1); yd=zeros(p,1); ye=zeros(p,1);
    xo=zeros(p,1); xd=zeros(p,1); xe=zeros(p,1);

    % get window params for default source distribution (x)
    f=floor(n/p); c=ceil(n/p); r=mod(n,p);

    % get window sizes and origins for default source distribution (x)
    od(1:p)=c;
    od(r+1:p)=f;
    oo(1)=1;
    for i=1:p-1; oo(i+1)=oo(i)+od(i); end;
    oe=oo+od-1;

    % get window params for target (y) oriented distribution
    m=n+(p-1)*2*h;
    f=floor(m/p); c=ceil(m/p); r=mod(m,p);

    % get window sizes and origins for target (y)
    yd(1:p)=c;
    yd(r+1:p)=f;
    yo(1)=1;
    for i=1:p-1; yo(i+1)=yo(i)+yd(i); end;
    ye=yo+yd-1;

    % get window sizes and origins for source (x)
    xd(1)=yd(1)-h;
    for i=2:p-1; xd(i)=yd(i)-2*h; end
    xd(p)=yd(p)-h;
    xo(1)=1;
    for i=1:p-1; xo(i+1)=xo(i)+xd(i); end;
    xe=xo+xd-1;

    % test for minimal window size for source
    t=min(xd);
    assert(t>2*h,'funWindow1DShape: half-overlap (%d) too large for local window size (%d). Args: (%d,%d,%d)',h,t,n,p,h);

    % put array together
    os = [oo od oe];
    ys = [yo yd ye];
    xs = [xo xd xe];

    % check if everything sums up
    cn=sum(xs(:,2));
    assert(cn==n,'funWindow1DShape: window sum %d != n=%d',cn,n);
    cm=sum(ys(:,2));
    assert(cm==m,'funWindow1DShape: extended window sum %d != m=%d',cm,m);

    % debuging only
    %disp('funWindow1DShape');
    %disp([n sum(cn) p h m sum(cm) c r f]);
    %disp([ys']); disp([xs']);

end
