function str = iwave_MPIcommand(nlab,flagReg);
if flagReg == 1
        str = ['mpiexec.hydra -np ' num2str(nlab)];
else
        str = ['mpiexec.hydra -np ' num2str(nlab) ' -hosts localhost'];
end

%str = ['mpirun -np ' num2str(nlab) ];

% ['mpiexec.hydra -np ' num2str(nlab) ' -hostfile ' getenv('PBS_NODEFILE')]
% ['mpirun -np num2str(nlab)'];

