function ps = parpool_size()
% parpool_size - return number of workers in parallel pool
% parpool_size returns 0 id parallel pool is not opened

dcv=ver('distcomp');
if str2double(dcv.Version) < 6.3
    ps = matlabpool('size');
else
    pool = gcp('nocreate');
    if isempty(pool)
        ps = 0;
    else
        ps = pool.NumWorkers;
    end
end

end
