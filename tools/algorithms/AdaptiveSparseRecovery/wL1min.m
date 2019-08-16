function [Dest] = wL1min(b, mask, options)
% Recover missing traces using weighted L1 minimization and utilizing the
% correlation between the support sets of the transform coefficients across
% partitions.
% 
%     b         - subsampled data volume
%     mask      - logical data volume indicating the locations of the missing data entries
%     options:
%         partOrder - specifies the proper permutation of the data volume so that 
%                    partitioning is always performed along the first dimension
%         st        - index of the first partition
%         fin       - index of the last partition
%         C	        - sparsifying transform
%         maxiter   - maximum number of spgl1 iterations
%         omega	    - sets the weight used in weighted L1 minimization

% Copyright 2013 Hassan Mansour (hassanm@cs.ubc.ca)


partOrder = getOption(options, 'partorder', [1 2 3]); % define the partition ordering
st        = getOption(options, 'st', 1);        % index of first partition
fin       = getOption(options, 'fin', size(b,1));     % index of last partition
C         = getOption(options, 'transform', 1); % sparsifying transform
maxiter   = getOption(options, 'maxiter', 500); % maximum number of spql1 iterations
omega     = getOption(options, 'omega', 0.3);   % initialize the weights


% permute the measurements so that partitioning is applied only across the
% first dimension
b    = permute(b,partOrder);
mask = permute(mask,partOrder);

dim = size(b);


% start with the zero offset, i.e. off = dim(2) and move to 1
x_prev=[];

opts.verbosity=1;
opts.optTol = 1e-4;
opts.bpTol = 1e-4;
opts.iterations=maxiter;
opts.weights = [];
        

x_prev=[];

if fin > st
    incr = 1;
else
    incr = -1;
end

for i=st:incr:fin
    if isempty(x_prev)
        
        n = size(C,2);
        N = size(C,1);
        
        y = vec(b(i,:));
        
        % build mask operator
        RM = opMask(dim(2)*dim(3),find(vec(mask(i,:)) > 0));
        
        % measurement matrix
        A = RM*C';
    
        sigma = 1e-4;
        
        W = ones(N,1);
        
        [x_l1] = spgl1(A, y, 0, sigma, zeros(N,1), opts);
        
        Dest(i,:,:) = reshape(C'*x_l1, dim(2), dim(3));
        x_prev = C*C'*x_l1;
       
        
    else
        n = size(C,2);
        N = size(C,1);

        y = vec(b(i,:));
        
        % build mask operator
        RM = opMask(dim(2)*dim(3),find(vec(mask(i,:)) > 0));

        [Cx Idx] = sort(abs(x_prev), 'descend');
        ratCx = sqrt(cumsum(Cx.^2))/norm(Cx);
        k = find(ratCx(1:end) >= 0.9);
        if isempty(k)
            T = Idx(1);
        else
            T = Idx(1:k(1));
        end
        Tc = setdiff([1:N]',T);

        % weight vector
        W = ones(N,1);
        W(T) = omega;
        
        % measurement matrix
        A = RM*C';

        sigma = 1e-4;

        % run weighted L1 minimization
        opts.weights = W;
        [x_wl1] = spgl1(A, y, 0, sigma, zeros(N,1), opts);

        Dest(i,:,:) = reshape(C'*x_wl1, dim(2), dim(3));
        x_prev = C*C'*x_wl1;

    end
end




    

