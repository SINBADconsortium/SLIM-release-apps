function pool = parpool_open(nw)
% parpool_open - opens default parallel pool of given size
% parpool_open returns empty 'N/A' if matalbppol is called
% or parallel-pool object if parpool is called

    pool = parpool(nw);
end
