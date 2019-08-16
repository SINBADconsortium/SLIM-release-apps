function pool = parpool_open(nw)
% parpool_open - opens default parallel pool of given size
% parpool_open returns empty 'N/A' if matalbppol is called
% or parallel-pool object if parpool is called

dcv=ver('distcomp');
if str2double(dcv.Version) < 6.3
    matlabpool('open',nw);
    pool = 'N/A';
else
    pool = parpool(nw);
end

end
