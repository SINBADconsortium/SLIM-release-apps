function x = distSort_timeSlice2shot(x)
% Takes in a distributed 2-d prestack seismic data sorted in shot gathers
% (3D-array of dt=recv_trace, d2=shot, and d3=time-slice) and re-sort into
% time-slice gathers (3D-array d1=time, d2=recv_trace, and d3=shot)
%    
%    essentially performs the following operation for distributed 3d-arrays:
%       x = permute(x, [3 1 2])
%
% (works with MATLAB Parallel Computing Toolbox v4.2)

    dims = size(x);
    nt = dims(3);
    nr = dims(1);
    ns = dims(2);

    % probe the behaviour of default distribution scheme using a test vector
    test_vec = distributed.zeros(ns,1);

    spmd
        codistr = getCodistributor(test_vec);
        partition_shot = codistr.Partition;

        new_codistr = codistributor1d(2, partition_shot, [nr ns nt]);
        new_codistr_permuted = codistributor1d(3, partition_shot, [nt nr ns]);

        % redistribute along first dim, then permute locally, then build as a 3d array at the end
        x = redistribute(x, new_codistr);
        local_x = getLocalPart(x);
        local_x = permute(local_x, [3 1 2]);
        x = codistributed.build(local_x, new_codistr_permuted,'noCommunication');
    end

end
