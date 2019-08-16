function op = opFFTsym_datacube(nt, ns, nr, mask, compress_size)
% OPFFTSYM_datacube  One-dimensional symmetric real fast Fourier transform (FFT).
%                        NOTE: Specialized version adapted for EPSI_SLIM
%
%    opFFTsym_datacube(nt, ns, nr) creates a one-dimensional
%    normalized Fourier transform along the first axis of Greens
%    function datacube
%
%    opFFTsym_datacube(nt, nr nr, mask) if mask is a 1D vector
%    creates a one-dimensional normalized Fourier transform along the
%    first axis of Greens function datacube, followed by a restriction
%    matrix in frequency defined by mask as an array of indicies
%    (note: need length(mask) <= nf, and max(mask) <= nf)
%
%    opFFTsym_datacube(nt, ns, nr, mask) if mask is a 2D vector
%    creates a one-dimensional normalized Fourier transform along the
%    first axis of Greens function datacube, followed by a
%    shot-dependent restriction matrix in frequency defined by mask as
%    an array selected frequency indecies on the first axis and the
%    shot positions on the second axis.
%
%    opFFTsym_datacube(nt, ns, nr, mask, compress_size) if the
%    flag compress_size is set to 1, then the operator will compress
%    the matrix to the size of the restricted data. If the flag is
%    zero, the operatro simply zeroes out the masked bits.  The
%    default option is 1.
%    
%
%   Copyright 2008, Tim Lin


% calculate operator size
if mod(nt,2) == 1
    nf = ceil(nt/2);
    has_nyquist = 0;
else
    nf = (nt/2) + 1;
    has_nyquist = 1;
end

m = nf * nr * ns;
n = nt * nr * ns;

if ~exist('mask','var')
    % no masking, normal FFT on datacube 1st dimension
    
    subfunc_handle = @(x,mode) opFFTsym_datacube_intrnl(nt,ns,nr,nf,m,mft,n,has_nyquist,x,mode);
    op = opFunction(m,n,subfunc_handle); % return a SPOT operator using constructor opFunction
    
else
    % Determine number of non-masked frequencies
    if islogical(mask)
        nft = sum(mask);
    else
        nft = size(mask,1);
    end
    
    % Determine whether to delete the masked frequencies
    if ~exist('compress_size','var')
        compress_size = 1;
    end
    
    if compress_size == 1
        if islogical(mask)
            mft = sum(mask)*ns*nr;
        else
            mft = size(mask,1)*ns*nr;
        end
    else
        mft = m;
    end
    
    subfunc_handle = @(x,mode) opFFTsym_datacube_intrnl(nt,ns,nr,nf,m,mft,n,has_nyquist,x,mode,mask,compress_size);
    op = opFunction(mft,n,subfunc_handle); % return a SPOT operator using constructor opFunction
end


function y = opFFTsym_datacube_intrnl(nt,ns,nr,nf,m,mft,n,has_nyquist,x,mode,mask,compress_size)

    % apply operator
    if mode == 0

        if ~exist('mask','var')
             % no mask defined
             y = {m,n,[0,0,0,0],{'FFTsym_datacube'}};
        else
             y = {mft,n,[0,0,0,0],{'FFTsym_datacube with frequency restrictions'}};
        end
   
    elseif mode == 1
       
       % reshape and apply fft
       x = reshape(x, nt, ns, nr);
       x = fft(x);
       y = x;
       clear x;
       y(nf+1:end,:,:) = [];
       y = y ./ sqrt(nt);
       % doubling the magnitude of everything except DC and Nyquist, because the conjugate freqs are chopped off
       y(2:end-1,:,:) = 2.*y(2:end-1,:,:);
       if not(has_nyquist), y(end,:,:) = 2.*y(end,:,:); end
       
       
       % check for mask and mask type
       if exist('mask') == 1 
            lm = length(mask(:));
            if lm <= nf
                % uniform frequency masking across the shots
                if compress_size == 1
                    y = y(mask,:,:);
                else
                    if islogical(mask)
                        y(~mask,:,:)  = 0;  % mask is stencil
                    else
                        y(setdiff([1:nf],mask),:,:) = 0; % mask is list of indices
                    end
                end
            elseif lm <= ns * nf
                % independent frequency masking across the shots
                if compress_size == 1
                    ytemp = zeros(size(mask,1),ns,nr);
                    for ks = 1:ns
                        ytemp(:,:,ks) = y((mask(:,ks)),:,ks);
                    end
                    y = ytemp;
                else
                    for ks = 1:ns
                        y(setdiff([1:nf],mask(:,ks)),:,ks) = 0;
                    end
                end
            else
                error('error: mask size incorrect');
            end
       end
       y = y(:);
       
    else % mode 2
    
    if exist('mask','var')
        if compress_size == 1
            % zero-fill the masked entries
            xtemp = zeros(nf,ns,nr,class(x));
            x = reshape(x,nft,ns,nr);

            % check whether restriction is uniform or shot-dependent
            if size(mask,2) > 1
                for ks = 1:ns
                    xtemp(mask(:,ks),:,ks) = x(:,:,ks);
                end
            else
                xtemp(mask,:,:) = x;
            end
            x = xtemp;
            clear xtemp   
        else
        x = reshape(x, nf, ns, nr);
        end
    else
        x = reshape(x, nf, ns, nr);
    end

       % symmetric flag automatically fills in the conjugate frequencies
       x = ifft(x,nt,'symmetric') .* sqrt(nt);
       x = x(:);
       y = x;
       clear x;
    end
end
end