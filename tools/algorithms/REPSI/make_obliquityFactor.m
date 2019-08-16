function [obfac invobfac] = make_obliquityFactor(dt,nf,frequency_mask)

    % make a list of omega
    f_max = 0.5 / dt;
    df = f_max / (nf - 1);
    om = 2*pi*[0:df:f_max];

    %build obliquity factor
    % Define obliquity factor in frequency domain
    obfac = (sqrt(om) * exp(sqrt(-1)*pi*0.25));
    invobfac = 1 ./ obfac;
    invobfac(1) = 0;
    
    if exist('frequency_mask','var')
        frequency_mask = logical(frequency_mask);
        obfac = obfac(frequency_mask);
        invobfac = invobfac(frequency_mask);
    end

    obfac = obfac(:);
    invobfac = invobfac(:);
    
end