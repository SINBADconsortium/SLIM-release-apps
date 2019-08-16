function parpool_close()
% parpool_close - closes parallel pool

dcv=ver('distcomp');
if str2double(dcv.Version) < 6.3
    matlabpool('close');
else
    delete(gcp('nocreate'));
end

end
