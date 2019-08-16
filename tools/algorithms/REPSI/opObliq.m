function op = opObliq(dt,nf_total,nd,frequency_mask)
% OPOBLIQ       appends obliquity factor in the Fourier domain
% dt is time interval (in s)
% nf_total is number of samples in frequency axis
% nd is number of traces in other directions


if exist('frequency_mask','var')
    [obfac invobfac] = make_obliquityFactor(dt,nf_total,frequency_mask);
else
    [obfac invobfac] = make_obliquityFactor(dt,nf_total);
end

subfunc_handle = @(x,mode) opObliq_intrnl(x,mode);
% op = @(x,mode) opObliq_intrnl(x,mode);
m = length(obfac)*nd;
n = length(obfac)*nd;

op = opFunction(m,n,subfunc_handle); % return a SPOT operator using constructor opFunction

    function y = opObliq_intrnl(x,mode);

    switch mode
        case 0
            y = {nf*nd,nf*nd,[0,0,0,0],{'opObliq'}};

        case 1
            sig = repmat(obfac,nd,1);
            x = sig.*x;
            y = x;
            clear sig

        case 2
            sig = repmat(conj(obfac),nd,1);
            x = sig.*x;
            y = x;
            clear sig
    end
    end
end
