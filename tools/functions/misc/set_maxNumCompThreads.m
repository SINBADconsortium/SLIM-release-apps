function [n N] = set_maxNumCompThreads(varargin)
% [n N] = set_maxNumCompThreads	controls the number of computational threads
%    on workers and master, i.e. maxNumCompThreads, and returns those values.
%
%    set_maxNumCompThreads sets values to MATLAB's defaults.
%
%    set_maxNumCompThreads(n) sets the number of computational threads
%    on workers and master to n.
%
%    set_maxNumCompThreads(n,N) sets the number of computational threads
%    on workers to n and on master to N.
%
%    set_maxNumCompThreads(n,'automatic') sets the number of computational threads
%    on workers to n and on master to number of cores.
%
%    NOTE: it is a BAD idea to set n to 'automatic' on workers.
%
    assert(nargin<=2,'Two optional argument at most')

    n=1;
    N='automatic';

    if nargin==1
        n=varargin{1};
        N=n;
    elseif nargin==2
        n=varargin{1};
        N=varargin{2};
    end
    
    spmd
        warning('off','MATLAB:maxNumCompThreads:Deprecated')
        maxNumCompThreads(n);
        setenv('OMP_NUM_THREADS',int2str(n))
    end
    spmd
        assert(maxNumCompThreads==n, 'Does not work')
	assert(str2num(getenv('OMP_NUM_THREADS'))==n, 'Does not work')
    end
    warning('off','MATLAB:maxNumCompThreads:Deprecated')
    maxNumCompThreads(N);
    N = maxNumCompThreads;
    setenv('OMP_NUM_THREADS',int2str(N))

end
