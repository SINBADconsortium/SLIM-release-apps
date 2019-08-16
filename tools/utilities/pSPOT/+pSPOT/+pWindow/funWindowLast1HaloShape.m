function [ m xs ys ] = funWindowLast1HaloShape( n, p, h )
%funWindowLast1HaloShape is a support function for forward oplWindow1D* operators.
%
%   [ M XS YS ] = funWindowLast1HaloShape( N, P, H )
%
%   INPUT:
%      N = length of the input vector
%      P = number of processors
%      H = half of the overlap's size
%   OUTPUT:
%      M = length of the output vector
%      XS = (p,3) vector holding start, size, end indecies
%           of the default distribution of the input vector in every window
%      YS = (p,3) vector holding start, size, end indecies
%           of the default distribution of the output vector in every window

    % check # of processors
    assert(p>1,'funWindowLast1HaloShape: number of processors has to be bigger then %d',p);

    % initialize shape vecs
    xo=zeros(p,1); xd=zeros(p,1); xe=zeros(p,1);
    yo=zeros(p,1); yd=zeros(p,1); ye=zeros(p,1);

    % get window params for default source distribution (x)
    f=floor(n/p); c=ceil(n/p); r=mod(n,p);

    % get window sizes and origins for default source distribution (x)
    xd(1:p)=c;
    xd(r+1:p)=f;
    xo(1)=1;
    for i=1:p-1; xo(i+1)=xo(i)+xd(i); end;
    xe=xo+xd-1;

    % get window params for target (y) oriented distribution
    m=n+(p-1)*2*h;

    % get window sizes and origins for target (y)
    yd(1:p)=c;
    yd(r+1:p)=f;
    yd(1)=yd(1)+h;
    yd(2:p-1)=yd(2:p-1)+2*h;
    yd(p)=yd(p)+h;
    yo(1)=1;
    for i=1:p-1; yo(i+1)=yo(i)+yd(i); end;
    ye=yo+yd-1;

    % test for minimal window size for source
    t=min(xd);
    assert(t>2*h,'funWindowLast1HaloShape: half-overlap (%d) too large for local window size (%d). Args: (%d,%d,%d)',h,t,n,p,h);

    % put array together
    ys = [yo yd ye];
    xs = [xo xd xe];

    % check if everything sums up
    cn=sum(xs(:,2));
    assert(cn==n,'funWindowLast1HaloShape: window sum %d != n=%d',cn,n);
    cm=sum(ys(:,2));
    assert(cm==m,'funWindowLast1HaloShape: extended window sum %d != m=%d',cm,m);

    % debuging only
    %disp('funWindowLast1HaloShape');
    %disp([n sum(cn) p h m sum(cm) c r f]);
    %disp([ys']); disp([xs']);

end
