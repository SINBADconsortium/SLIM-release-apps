function [D, J] = Fd2DT(m,model,options)
% Time domain FD modeling operator by iWAVE
%
% USE:
%    [D, J] = Fd2DT(m,model,options)
% input:
%   m                 - vector with gridded parameters [bulk [10^9 kg/(m*s^2)]; buoyancy [cm^3/gr]]
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc} - vectors describing source array
%   model.{zrec,xrec} - vectors describing receiver array.
%   options           - options for iWAVE
%   options.fwdpara   - the forward operator parameter file
%   options.adjpara   - the adjoint operator parameter file
%   options.linpara   - the linear operator parameter file
%   options.delete      - whether delete or not delete the file that iWave creates
%
% output:
%   D  - Data cube (nrec x nsrc x ntime)
%   J  - Jacobian as an operator

% Author: Zhilong Fang
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: August, 2013
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if isfield(model,'randseq')
    isrand = model.randseq;
else
    isrand = 0;
end

if isfield(model,'simultsrc')
    ismulti = model.simultsrc;
else
    ismulti = 0;
end

curdir = pwd;
%addpath(curdir);
if isfield(options,'label')
    label = [options.label,'_fwd'];
else
    label = 'fwd_results_temp';
end
model.label = label;
mkdir(label);
% Read bulk modulus and buoyancy
bulk = reshape(m(1:length(m)/2),model.n);
buoy = reshape(m(1+length(m)/2:end),model.n);

% Write bulk modules and buoyancy into an rsf file
d   = model.d;
o   = model.o;
n   = model.n;
% WriteAllData([curdir '/' label '_buoy.rsf'],buoy,n,d,o);
% WriteAllData([curdir '/' label '_bulk.rsf'],bulk,n,d,o);
odnwrite([curdir '/' label '_buoy.rsf'],buoy,o,d,n);
odnwrite([curdir '/' label '_bulk.rsf'],bulk,o,d,n);

fname_buoy = [curdir '/' label '_buoy.rsf'];
fname_bulk = [curdir '/' label '_bulk.rsf'];

% Modify the rsf header file
fid = fopen(fname_buoy,'a');
str = ['data_type=buoyancy\n'];
fprintf(fid,str);
%str = ['in=' curdir '/' label '_buoy.rsf@\n'];
%fprintf(fid,str);
fclose(fid);

fid = fopen(fname_bulk,'a');
str = ['data_type=bulkmod\n'];
fprintf(fid,str);
%str = ['in=' curdir '/' label '_bulk.rsf@\n'];
%fprintf(fid,str);
fclose(fid);

fwd_params = importdata(options.fwdpara);
fname_buoy_temp = ['../' label '_buoy.rsf'];
fname_bulk_temp = ['../' label '_bulk.rsf'];

if isrand == 0 && ismulti == 0
    %% Forward part
    fwd_params{end+1} = 'srctype=''point'';';
    [status,result]=fwdwrapper(fname_buoy_temp,fname_bulk_temp,model,fwd_params,options);
    assert(status==0,'iwave++ fails');
    data = [label '_data.su'];
    nt   = length(model.t);
    D    = ReadSuFast(data,nt,'','b');
    
    cd(curdir);
    
    if isfield(options,'rand')
        if isfield(options,'delete')
            if options.delete == 1
                rmf            = ['rm -rf ' label '*'];
                [status,result]=system(rmf);
                delete(['cout*']);
            end
        end
        
    else
        if isfield(options,'delete')
            if options.delete == 1
                str = ['!rm -rf ' label];
                eval(str);
                %rmdir(label,'s');
                rmf            = ['rm -rf ' label '*'];
                [status,result]=system(rmf);
                %        			delete([label '*']);
            end
        end
    end
else
    % Modification for parallel computing other the sources(vts)
    if isrand > 0
        Numsrc = length(model.xsrc);
    else
        Numsrc = size(model.srcdata,3);
    end
        
    
    
    spmd
        
        workidx         = labindex;
        modellocal      = model;
        Part            = codistributor1d.defaultPartition(Numsrc);
        Numsrc_local    = Part(workidx);
        nlab            = numlabs;
        
        fname_buoy_temp = [ label '_buoy.rsf'];
        fname_bulk_temp = [ label '_bulk.rsf'];
        
        % Repartition of the sources between the workers
        [label_temp]             = [label '_local' num2str(workidx)];
        id_localsrc              = sum(Part(1:workidx-1))+1:sum(Part(1:workidx));
        if isrand > 0
            modellocal.xsrc          = modellocal.xsrc(id_localsrc);
        else
            modellocal.srcdata       = modellocal.srcdata(:,:,id_localsrc);
        end
        fwd_params_locali        = importdata(options.fwdpara);
        fwd_params_locali{end+1} = 'partask=1;';
        if isrand == 0
            fwd_params_locali{end+1} = 'srctype=''array'';';
        end
        D_local                  = [];
        
        for i=1:Numsrc_local
            
            modellocali             = modellocal;
            if isrand > 0
                modellocali.xsrc    = modellocal.xsrc(i);
            else
                modellocali.srcdata = squeeze(modellocal.srcdata(:,:,i));
            end
            options_locali          = options;
            label_locali            = [label_temp '_' num2str(i)]; % The adjpara.txt should be modified
            options_locali.label    = label_locali;
            
            fname_buoy_temp_locali  = [label_locali '_buoy.rsf'];
            fname_bulk_temp_locali  = [label_locali '_bulk.rsf'];
            
            
            [status,result]         = fwdwrapper(fname_buoy_temp,fname_bulk_temp,modellocali,fwd_params_locali,options_locali);
            assert(status==0,'iwave++ fails');
            data_label_locali       = [label_locali '_fwd_data.su'];
            D_locali                = ReadSuFast(data_label_locali,length(model.t),'','b');
            D_local                 = [D_local;D_locali(:)];
        end
        
        Part   = codistributor1d.defaultPartition(Numsrc);
        Part   = Part * length(model.t) * length(model.xrec);
        codist = codistributor1d(1, Part, [sum(Part),1]);
        D      = codistributed.build(D_local,codist,'noCommunication');
    end %end spmd
    
    D = gather(D);
    D = reshape(D,length(model.t),length(D)/length(model.t));
    
    if isfield(options,'delete')
        if options.delete == 1
            rmf            = ['rm -rf ' label '*'];
            [status,result]=system(rmf);
            delete(['cout*']);
        end
    end
end

% Construct SPOT operator
if nargout > 1
    J = oppDFd2DT(m,model,options);
end
end




