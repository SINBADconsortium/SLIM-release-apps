function parpool_close()
% parpool_close - closes parallel pool

    delete(gcp('nocreate'));
end
