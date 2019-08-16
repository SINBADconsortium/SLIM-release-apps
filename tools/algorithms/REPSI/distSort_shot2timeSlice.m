function x = distSort_shot2timeSlice(x)
% Takes in a distributed 2-d prestack seismic data sorted in shot gathers
% (3D-array d1=time, d2=recv_trace, and d3=shot) and re-sort into
% time-slice gathers (3D-array of dt=recv_trace, d2=shot, and d3=time-slice)
%    
%    essentially performs the following operation for distributed 3d-arrays:
%       x = permute(x, [2 3 1])
%
% (works with MATLAB Parallel Computing Toolbox v4.2)

    dims = size(x);
    nt = dims(1);
    nr = dims(2);
    ns = dims(3);

    % probe the behaviour of default distribution scheme using a test vector
    test_vec = distributed.zeros(nt,1);
    
    spmd
        codistr = getCodistributor(test_vec);
        partition_time = codistr.Partition;
        
        new_codistr = codistributor1d(1, partition_time, [nt nr ns]);
        new_codistr_permuted = codistributor1d(3, partition_time, [nr ns nt]);
        
        % redistribute along first dim, then permute locally, then build as a 3d array at the end
        x = redistribute(x, new_codistr);
        local_x = getLocalPart(x);
        local_x = permute(local_x, [2 3 1]);
        x = codistributed.build(local_x, new_codistr_permuted,'noCommunication');
    end
    
end
